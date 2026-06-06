#!/usr/bin/env python3
"""
asd_duplex.py

Workflow:
1. Scan all k-mer windows and score against their exact reverse complement.
2. Report ALL raw windows.
3. Find all windows tied for the best ΔG (within tolerance).
4. Merge overlapping/adjacent top-hit windows.
5. Extract merged sequences from the source RNA.
6. Re-score merged sequences.
7. Report merged candidates separately.
"""

from dataclasses import dataclass
import argparse
import RNA

DEFAULT_TEMP = 37.0
DEFAULT_WINDOW_K = 4
DEFAULT_ASD_TAIL_NT = 25
DG_TOL = 1e-6


def load_fasta_sequence(path, to_rna=True):
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

# returns the complementary RNA seq
def revcomp_rna(seq):
    comp = str.maketrans("AUCG", "UAGC")
    return seq.upper().replace("T", "U").translate(comp)[::-1]


def duplex_dg(seq1, seq2, temp=DEFAULT_TEMP):
    RNA.cvar.temperature = temp
    return RNA.duplexfold(seq1, seq2).energy


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


def window_scan(seq, abs_offset, k, temp):
    hits = []
    for i in range(len(seq) - k + 1):
        sub = seq[i:i+k]
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


def merge_top_hits(hits, source_seq, source_start, temp, tolerance=DG_TOL):
    if not hits:
        return []

    best_dg = hits[0].dg

    top_hits = [
        h for h in hits
        if abs(h.dg - best_dg) <= tolerance
    ]

    top_hits.sort(key=lambda h: h.abs_start)

    merged_intervals = []
    cur_start = top_hits[0].abs_start
    cur_end = top_hits[0].abs_end

    for h in top_hits[1:]:
        if h.abs_start <= cur_end + 1:
            cur_end = max(cur_end, h.abs_end)
        else:
            merged_intervals.append((cur_start, cur_end))
            cur_start, cur_end = h.abs_start, h.abs_end

    merged_intervals.append((cur_start, cur_end))

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


def build_parser():
    p = argparse.ArgumentParser()

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence")
    #src.add_argument("--asd-seq") #add seq will only work if func are modified not to take end/start/tail but also accept free seq only

    p.add_argument("--asd-start", type=int)
    p.add_argument("--asd-end", type=int)
    p.add_argument("--asd-tail-nt", type=int, default=DEFAULT_ASD_TAIL_NT)

    p.add_argument("--window-k", type=int, default=DEFAULT_WINDOW_K)
    p.add_argument("--temp", type=float, default=DEFAULT_TEMP)

    return p


def main():
    args = build_parser().parse_args()

    if args.asd_seq:
        seq = args.asd_seq.upper().replace("T", "U")
        abs_start = 1

    else:
        _, full_seq = load_fasta_sequence(args.sequence, True)

        end = args.asd_end or len(full_seq)
        start = args.asd_start or (len(full_seq) - args.asd_tail_nt + 1)

        seq = full_seq[start - 1:end]
        abs_start = start

    hits = window_scan(
        seq,
        abs_start,
        args.window_k,
        args.temp
    )

    print("RAW WINDOW HITS")
    print("ABS_START\tABS_END\tSUBSEQ\tDG")

    for h in hits:
        print(f"{h.abs_start}\t{h.abs_end}\t{h.subseq}\t{h.dg:.2f}")

    merged = merge_top_hits(
        hits,
        seq,
        abs_start,
        args.temp
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
