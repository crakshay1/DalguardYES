#!/usr/bin/env python3
"""
asd_duplex.py

Purpose:
    Identifies candidate anti-Shine-Dalgarno (ASD) sequences within an RNA
    molecule by scanning k-mer windows and scoring how well each window
    base-pairs with its own reverse complement (i.e. how good a self-
    complementary duplex it can form).  

    This script finds the region of the input RNA most likely to act as an ASD by
    looking for the segment with the most favourable (most negative) folding
    free energy (ΔG) when paired with its reverse complement.

Workflow:
    1. Scan all k-mer windows and score each against its exact reverse
       complement using ViennaRNA's duplexfold().
    2. Report ALL raw window scores (tab-delimited to stdout).
    3. Find all windows tied for the best ΔG (within a floating-point
       tolerance), i.e. the globally optimal windows.
    4. Merge overlapping or adjacent top-hit windows into contiguous regions.
    5. Extract those merged regions from the original RNA sequence.
    6. Re-score each merged region as a whole duplex.
    7. Report merged candidates separately.

@cuajiniquil (EmmaGH)
"""

from dataclasses import dataclass  
import argparse                     
import RNA                          

# ─────────────────────────────────────────────────────────────────────────────
# Module-level constants / defaults

DEFAULT_TEMP = 37.0        # Folding temperature in °C
DEFAULT_WINDOW_K = 4       # Default k-mer length for the sliding window scan
DEFAULT_ASD_TAIL_NT = 25   # How many nucleotides from the 3' end to analyse when no --asd-start/--asd-end is given.
                           # 25 nt covers the typical ASD region in 16S rRNA.
DG_TOL = 1e-6              
# Floating-point tolerance for comparing ΔG values. 
# Two energies within 1×10⁻⁶ kcal/mol are treated as identical (avoids fp rounding artefacts).


# ─────────────────────────────────────────────────────────────────────────────
# I/O helper

def load_fasta_sequence(path, to_rna=True):
    """
        Read the first sequence record from a FASTA file.
        Returns the full concatenated sequence, upper-cased with its header.
    """
    header, parts = "", []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue                    # skip blank lines
            if line.startswith(">"):
                if parts:
                    break                   # second record found → stop
                header = line[1:]           # strip the leading '>'
            else:
                seq = line.upper()          # normalise to upper-case
                if to_rna:
                    seq = seq.replace("T", "U")   # DNA → RNA
                parts.append(seq)
    return header, "".join(parts)


# ─────────────────────────────────────────────────────────────────────────────
# Sequence utilities

def revcomp_rna(seq):
    """
        Return the reverse complement of an RNA sequence.
        Watson-Crick RNA base-pairing:  A↔U,  U↔A,  C↔G,  G↔C  

        1. Upper-case and convert any 'T' → 'U' (tolerates accidental DNA input).
        2. Translate each base to its complement using str.maketrans / translate.
        3. Reverse the resulting complement string.

        Example :
        revcomp_rna("AUGC") → complement "UACG" → reversed "GCAU"
    """
    comp = str.maketrans("AUCG", "UAGC")          # maps A→U, U→A, C→G, G→C
    return seq.upper().replace("T", "U").translate(comp)[::-1]


# ─────────────────────────────────────────────────────────────────────────────
# Thermodynamic scoring

def duplex_dg(seq1, seq2, temp=DEFAULT_TEMP):
    """
        Compute the hybridisation free energy (ΔG, kcal/mol) of an RNA duplex
        formed by seq1 and seq2 using ViennaRNA's duplexfold algorithm.  
        A more negative ΔG indicates a more stable (stronger) duplex.  
        Returns ΔG in kcal/mol.
    """
    RNA.cvar.temperature = temp           # set ViennaRNA's global temperature
    return RNA.duplexfold(seq1, seq2).energy


# ─────────────────────────────────────────────────────────────────────────────
# Result data containers

@dataclass
class WindowHit:
    """
        Stores the result of scoring a single k-mer window.
    """

    abs_start: int  # 1-based start position in the full input sequence.
    abs_end: int    # 1-based end position (inclusive).
    subseq: str     # The k-mer RNA sequence.
    dg: float       # Duplex ΔG (kcal/mol) of subseq paired with its own reverse complement.


@dataclass
class MergedHit:
    """
        Stores the result of re-scoring a merged (multi-window) region.
    """
    abs_start: int  # 1-based start position in the full input sequence.
    abs_end: int    # 1-based end position (inclusive).
    subseq: str     # The full merged RNA sequence.
    dg: float       # Duplex ΔG of the entire merged region paired with its own reverse complement.


# ─────────────────────────────────────────────────────────────────────────────
# Core algorithm – sliding window scan

def window_scan(seq, abs_offset, k, temp):
    """
        Slide a window of length k across `seq` and score each k-mer.
        For each window position i the k-mer is scored against its own reverse
        complement:  ΔG = duplex_dg(seq[i:i+k], revcomp_rna(seq[i:i+k]))
        Returns List[WindowHit] – All windows, sorted ascending by ΔG (best first).
    """
    hits = []
    for i in range(len(seq) - k + 1):      # iterate over all valid windows
        sub = seq[i:i+k]
        hits.append(
            WindowHit(
                abs_start = abs_offset + i,         # 1-based absolute start
                abs_end   = abs_offset + i + k - 1, # 1-based absolute end
                subseq    = sub,
                dg        = duplex_dg(sub, revcomp_rna(sub), temp),
            )
        )
    hits.sort(key=lambda h: h.dg)   # sort: most stable (lowest ΔG) first
    return hits


# ─────────────────────────────────────────────────────────────────────────────
# Core algorithm – merge top hits into contiguous regions

def merge_top_hits(hits, source_seq, source_start, temp, tolerance=DG_TOL):
    """
        Collect all windows tied for the best ΔG, merge any that are
        adjacent or overlapping into contiguous intervals, then re-score each
        merged interval.

        1. Identify the best ΔG value (hits[0].dg, since `hits` is pre-sorted).
        2. Collect all hits whose ΔG differs by at most `tolerance` from best_dg.
        3. Sort those top hits by abs_start.
        4. Sweep through them with a greedy interval-merge:
            - Maintain a "current interval" [cur_start, cur_end].
            - If the next hit's start ≤ cur_end + 1, extend cur_end.
            - Otherwise, close the current interval and start a new one.
        5. For each merged interval, extract the subsequence from `source_seq`,
        compute a fresh duplex ΔG, and build a MergedHit.
        6. Return MergedHit list sorted by ΔG.

        Returns List[MergedHit] – Merged regions, sorted by ΔG (best first).
    """
    if not hits:
        return []

    best_dg = hits[0].dg     # lowest (most negative) energy in the sorted list

    # Collect all top-tier hits
    top_hits = [
        h for h in hits
        if abs(h.dg - best_dg) <= tolerance   # within floating-point tolerance
    ]

    # Sort by genomic position for interval sweep
    top_hits.sort(key=lambda h: h.abs_start)

    # Greedy interval merge
    merged_intervals = []
    cur_start = top_hits[0].abs_start
    cur_end   = top_hits[0].abs_end

    for h in top_hits[1:]:
        if h.abs_start <= cur_end + 1:
            # Overlapping or adjacent → extend the current interval
            cur_end = max(cur_end, h.abs_end)
        else:
            # Gap found → close current interval, start a fresh one
            merged_intervals.append((cur_start, cur_end))
            cur_start, cur_end = h.abs_start, h.abs_end

    merged_intervals.append((cur_start, cur_end))   # close the last interval

    # Extract sequences and re-score
    merged = []
    for start, end in merged_intervals:
        # Convert absolute 1-based coordinates to 0-based indices into source_seq
        rel_start = start - source_start       # offset from source_seq[0]
        rel_end   = end   - source_start + 1   # +1 for slice end

        seq = source_seq[rel_start:rel_end]

        merged.append(
            MergedHit(
                abs_start = start,
                abs_end   = end,
                subseq    = seq,
                # Re-score the entire merged region (may differ from the
                # individual k-mer scores because the region is longer)
                dg = duplex_dg(seq, revcomp_rna(seq), temp),
            )
        )

    merged.sort(key=lambda x: x.dg)   # best ΔG first
    return merged


# ─────────────────────────────────────────────────────────────────────────────
# CLI argument parser

def build_parser():
    p = argparse.ArgumentParser()
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence")
    p.add_argument("--asd-start",   type=int)
    p.add_argument("--asd-end",     type=int)
    p.add_argument("--asd-tail-nt", type=int, default=DEFAULT_ASD_TAIL_NT)
    p.add_argument("--window-k", type=int,   default=DEFAULT_WINDOW_K)
    p.add_argument("--temp",     type=float, default=DEFAULT_TEMP)

    return p


# ─────────────────────────────────────────────────────────────────────────────
# Execution

def main():
    args = build_parser().parse_args()
    if args.asd_seq:
        seq = args.asd_seq.upper().replace("T", "U")
        abs_start = 1   

    
    else:
        _, full_seq = load_fasta_sequence(args.sequence, True)
        end   = args.asd_end   or len(full_seq)
        start = args.asd_start or (len(full_seq) - args.asd_tail_nt + 1)

        seq       = full_seq[start - 1 : end]
        abs_start = start 

    hits = window_scan(
        seq,          
        abs_start,    
        args.window_k,
        args.temp,
    )
    print("RAW WINDOW HITS")
    print("ABS_START\tABS_END\tSUBSEQ\tDG")  

    for h in hits:
        print(f"{h.abs_start}\t{h.abs_end}\t{h.subseq}\t{h.dg:.2f}")

    merged = merge_top_hits(
        hits,
        seq,          
        abs_start,    
        args.temp,
    )

    print("MERGED TOP-HIT REGIONS")
    print("ABS_START\tABS_END\tLEN\tSEQUENCE\tDG")   

    for m in merged:
        print(
            f"{m.abs_start}\t{m.abs_end}\t"
            f"{len(m.subseq)}\t{m.subseq}\t{m.dg:.2f}"
        )


if __name__ == "__main__":
    main()