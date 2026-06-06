#!/usr/bin/env python3
"""
orbs_duplex.py

Scans an input RNA transcript, calculates the duplex energy for each potential self-complementary 
o-RBS candidate, print the scan results, and write the best merged* hit and core hit
to a standardized JSON file as <gene>_rbs_core.json.

*the best hit is the one with the lowest (most negative and most stable) duplex energy;
if the best (most negative) hit has adjacent hits with similar duplex energy 
(delta dg <= DG_TOL) are merged into a contiguous region

Uses: 

python3 orbs_duplex.py --sequence AUCGCAAAACAGGAUCGUC --mrna5 AAAAACAAAAA --mrna3 UUUUUUGUUUUUU
(output: query_20260606_201310_rbs_core.json)

python3 orbs_duplex.py --fasta Tv3test.fa --mrna5 AAAAACAAAAA --mrna3 UUUUUUGUUUUUU
(output: query_20260606_224031_rbs_core.json)
"""

from __future__ import annotations
from datetime import datetime

from dataclasses import dataclass
import argparse
import json
from pathlib import Path
import re

import RNA

DEFAULT_TEMP = 37.0
DEFAULT_WINDOW_K = 4
DEFAULT_ASD_TAIL_NT = 12
DG_TOL = 1e-6
DEFAULT_NAME = datetime.now().strftime("query_%Y%m%d_%H%M%S")

def normalize_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def sanitize_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return cleaned.strip("._-") or "sequence"

# gets the first FASTA record sequence
def load_fasta_sequence(path: str, to_rna: bool = True) -> tuple[str, str]:
    header, parts = "", []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if parts:
                    break
                header = line[1:]
            else:
                seq = line.upper()
                if to_rna:
                    seq = seq.replace("T", "U")
                parts.append(seq)
    return header, "".join(parts)


def write_fasta(path: Path, header: str, sequence: str, width: int = 80) -> None:
    sequence = normalize_rna(sequence)
    with path.open("w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i:i + width] + "\n")


def write_fasta_records(path: Path, records: list[tuple[str, str]], width: int = 80) -> None:
    with path.open("w") as fh:
        for header, sequence in records:
            sequence = normalize_rna(sequence)
            fh.write(f">{header}\n")
            for i in range(0, len(sequence), width):
                fh.write(sequence[i:i + width] + "\n")


def resolve_input_sequence(raw_sequence: str | None, fasta_path: str | None) -> tuple[str, str]:
    if fasta_path:
        header, sequence = load_fasta_sequence(fasta_path, True)
        base_name = sanitize_name(header.split()[0] if header else Path(fasta_path).stem)
        return sequence, base_name
    if raw_sequence:
        return normalize_rna(raw_sequence), ""
    raise ValueError("Provide either --fasta or --sequence.")


def revcomp_rna(seq: str) -> str:
    comp = str.maketrans("AUCG", "UAGC")
    return normalize_rna(seq).translate(comp)[::-1]

# RNA energy calculation for seq1 (o-RBS candidate) and seq2 (reverse complement of the candidate)
# uses ViennaRNA's duplexfold function, which calculates the minimum free energy of the duplex formed by two RNA sequences
def duplex_dg(seq1: str, seq2: str, temp: float = DEFAULT_TEMP) -> float:
    RNA.cvar.temperature = temp
    return RNA.duplexfold(normalize_rna(seq1), normalize_rna(seq2)).energy


# dataclasses for storing window scan hits and merged hits as objects 
# absolute coordinates being relative to the full input sequence
@dataclass
class WindowHit:
    abs_start: int
    abs_end: int
    subseq: str
    dg: float


@dataclass
class MergedHit:
    abs_start: int
    abs_end: int
    subseq: str
    dg: float

"""

slides a window of size k across the input sequence, calculates the duplex energy 
for the k-mer and its reverse complement for each window subsequence.
Returns a list of WindowHit objects sorted by increasing duplex energy 
(so the best/most stable hits are first)

"""
def window_scan(seq: str, abs_offset: int, k: int, temp: float) -> list[WindowHit]:
    hits = []
    for i in range(len(seq) - k + 1):
        sub = seq[i:i + k]
        hits.append(
            WindowHit(
                abs_start=abs_offset + i,
                abs_end=abs_offset + i + k - 1,
                subseq=sub,
                dg=duplex_dg(sub, revcomp_rna(sub), temp),
            )
        )
    hits.sort(key=lambda h: h.dg)
    return hits

"""

collects the top hits with similar duplex energy (delta dg <= DG_TOL), 
merges them if they're adjacent/overlapping into a contiguous region,
keeping and updating the interval coordinates and subsequence for the merged region) and returns a list of MergedHit objects sorted by increasing duplex energy
returns MergedHit objects with absolute coordinates relative to the full input sequence

"""
def merge_top_hits(hits: list[WindowHit], source_seq: str, source_start: int, temp: float, tolerance: float = DG_TOL) -> list[MergedHit]:
    if not hits:
        return []

    best_dg = hits[0].dg
    top_hits = [h for h in hits if abs(h.dg - best_dg) <= tolerance]
    # rearrange top hits by absolute start coordinate to help merging 
    top_hits.sort(key=lambda h: h.abs_start)

    #merge relevant top hits into contiguous regions, keeping track of their absolute coordinates 
    merged_intervals: list[tuple[int, int]] = []
    cur_start = top_hits[0].abs_start
    cur_end = top_hits[0].abs_end

    for h in top_hits[1:]:
        if h.abs_start <= cur_end + 1:
            cur_end = max(cur_end, h.abs_end)
        else: # find gap: close the interval and start a new one
            merged_intervals.append((cur_start, cur_end))
            cur_start, cur_end = h.abs_start, h.abs_end

    merged_intervals.append((cur_start, cur_end))

    # recalculate the duplex energy for the merged regions, using the full merged subsequence and its reverse complement, 
    # store as MergedHit objects with absolute coordinates
    merged = []
    for start, end in merged_intervals:
        rel_start = start - source_start
        rel_end = end - source_start + 1
        seq = source_seq[rel_start:rel_end]
        merged.append(
            MergedHit(
                abs_start=start,
                abs_end=end,
                subseq=seq,
                dg=duplex_dg(seq, revcomp_rna(seq), temp),
            )
        )

    merged.sort(key=lambda x: x.dg)
    return merged


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Identify and export the best o-RBS candidate as JSON.")

    # fixed/unused arguments (to be removed/replaced in future edits)
    parser.add_argument("--asd-start", type=int, help="1-based start of the ASD search region.")
    parser.add_argument("--asd-end", type=int, help="1-based end of the ASD search region.")
    parser.add_argument("--output-dir", default=".", help="Directory where the <gene>_rbs_core.json file is written.")

    # o-RBS input: exactly one of --fasta or --sequence (mutually exclusive)
    orbs_source = parser.add_mutually_exclusive_group(required=True)
    orbs_source.add_argument("--fasta", help="o-RBS transcript input as a FASTA file.")
    orbs_source.add_argument("--sequence", help="o-RBS transcript input as a raw RNA/DNA sequence.")

    # flanking mRNA fragments: raw sequence or FASTA, mutually exclusive per flank
    mrna5_source = parser.add_mutually_exclusive_group(required=True)
    mrna5_source.add_argument("--mrna5", help="5' mRNA fragment upstream of the RBS; raw sequence.")
    mrna5_source.add_argument("--mrna5-fasta", help="5' mRNA fragment upstream of the RBS; FASTA file.")

    mrna3_source = parser.add_mutually_exclusive_group(required=True)
    mrna3_source.add_argument("--mrna3", help="3' mRNA fragment downstream of the RBS (CDS start); raw sequence.")
    mrna3_source.add_argument("--mrna3-fasta", help="3' mRNA fragment downstream of the RBS (CDS start); FASTA file.")

    # canonical RBS input: raw sequence or FASTA, mutually exclusive
    crbs_source = parser.add_mutually_exclusive_group()  # not required
    crbs_source.add_argument("--canonical-rbs", default="AUUCCUCCACUAG", help="Canonical RBS sequence; raw sequence.")
    crbs_source.add_argument("--canonical-rbs-fasta", help="Canonical RBS sequence; FASTA file.")

    # customizable arguments (with defaults)
    parser.add_argument("--asd-tail-nt", type=int, default=DEFAULT_ASD_TAIL_NT)
    parser.add_argument("--window-k", type=int, default=DEFAULT_WINDOW_K)
    parser.add_argument("--temp", type=float, default=DEFAULT_TEMP)
    parser.add_argument("--name", default=DEFAULT_NAME, help="Gene name used for output naming.")

    return parser


def main() -> None:
    args = build_parser().parse_args()

    # resolver functions/if statements for fasta/raw sequence inputs for both RBS and flanking mRNA fragments, 
    # returning the sequence and an optional header (for name inference)
    full_seq, inferred_name = resolve_input_sequence(args.sequence, args.fasta)

    if args.mrna5_fasta:
        _, mrna5 = load_fasta_sequence(args.mrna5_fasta, to_rna=True)
    else:
        mrna5 = normalize_rna(args.mrna5)

    if args.mrna3_fasta:
        _, mrna3 = load_fasta_sequence(args.mrna3_fasta, to_rna=True)
    else:
        mrna3 = normalize_rna(args.mrna3)

    if args.canonical_rbs_fasta:
        _, canonical_rbs = load_fasta_sequence(args.canonical_rbs_fasta, to_rna=True)
    else:
        canonical_rbs = normalize_rna(args.canonical_rbs)

    gene_name = sanitize_name(args.name or inferred_name)

    end = args.asd_end or len(full_seq)
    start = args.asd_start or (len(full_seq) - args.asd_tail_nt + 1)
    if start < 1 or end > len(full_seq) or start > end:
        raise SystemExit("Invalid ASD window coordinates.")

    seq = full_seq[start - 1:end]
    abs_start = start

    hits = window_scan(seq, abs_start, args.window_k, args.temp)
    print("RAW WINDOW HITS")
    print("ABS_START\tABS_END\tSUBSEQ\tDG")
    for h in hits:
        print(f"{h.abs_start}\t{h.abs_end}\t{h.subseq}\t{h.dg:.2f}")

    merged = merge_top_hits(hits, seq, abs_start, args.temp)
    print("MERGED TOP-HIT REGIONS")
    print("ABS_START\tABS_END\tLEN\tSEQUENCE\tDG")

    for m in merged:
        print(f"{m.abs_start}\t{m.abs_end}\t{len(m.subseq)}\t{m.subseq}\t{m.dg:.2f}")

    # write output in format compatible with downstream scripts, 
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{gene_name}_rbs_core.json"
        best = merged[0]
        core = hits[0]
        candidate = {
            "name": gene_name,
            "five_prime_flank": mrna5,
            "canonical_rbs": canonical_rbs, # aka wt rbs, wild type rbs 
            "rbs": revcomp_rna(full_seq),
            "core": revcomp_rna(best.subseq),
            "spacer": "",
            "cds_start": mrna3,
            "mutable_regions": ["rbs", "five_prime_flank", "spacer"],
        }
        with output_path.open("w") as fh:
            json.dump([candidate], fh, indent=4)
            fh.write("\n")
        print(f"WROTE {output_path}")
    else:
        print("No merged o-RBS candidate found, no output written.")


if __name__ == "__main__":
    main()