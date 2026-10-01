#!/usr/bin/env python3

"""

orbs_duplex.py



Scans an input RNA transcript / anti-SD sequence for the strongest

self-complementary o-RBS core.



Important design rule used by the DalGuardYES pipeline:

- the strongest window hit defines the protected RBS core;

- the FULL o-RBS seed is the reverse complement of the entire ASD scan window,

  not only the strongest hit;

- therefore the core can have mutable RBS sequence on its left and right.



The output seed uses the project schema:

    five_prime_utr   fixed

    standby          mutable

    rbs_left         mutable

    rbs_core         fixed

    rbs_right        mutable

    spacer           mutable downstream

    cds_start        fixed



Example:

    python3 orbs_duplex.py \\

        --sequence ACTTGTATA \\

        --mrna5 CTACTAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA \\

        --standby AAAA \\

        --mrna3 ATGCGTAAAGGCGAAGAACTGTTTACCGGTGTGGTT \\

        --name p540-test-1

"""



from __future__ import annotations



import argparse

import json

import re

from dataclasses import dataclass

from datetime import datetime

from pathlib import Path



import RNA





DEFAULT_TEMP = 37.0

DEFAULT_WINDOW_K = 4

DEFAULT_ASD_TAIL_NT = 12

DG_TOL = 1e-6

DEFAULT_NAME = datetime.now().strftime("query_%Y%m%d_%H%M%S")





def normalize_rna(seq: str) -> str:

    """Uppercase a sequence and convert DNA T to RNA U."""

    return (seq or "").upper().replace("T", "U")





def sanitize_name(name: str) -> str:

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())

    return cleaned.strip("._-") or "sequence"





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





def write_fasta_records(

    path: Path,

    records: list[tuple[str, str]],

    width: int = 80,

) -> None:

    with path.open("w") as fh:

        for header, sequence in records:

            sequence = normalize_rna(sequence)

            fh.write(f">{header}\n")

            for i in range(0, len(sequence), width):

                fh.write(sequence[i:i + width] + "\n")





def resolve_input_sequence(

    raw_sequence: str | None,

    fasta_path: str | None,

) -> tuple[str, str]:

    if fasta_path:

        header, sequence = load_fasta_sequence(fasta_path, True)

        base_name = sanitize_name(

            header.split()[0] if header else Path(fasta_path).stem

        )

        return sequence, base_name



    if raw_sequence:

        return normalize_rna(raw_sequence), ""



    raise ValueError("Provide either --fasta or --sequence.")





def revcomp_rna(seq: str) -> str:

    comp = str.maketrans("AUCG", "UAGC")

    return normalize_rna(seq).translate(comp)[::-1]





def duplex_dg(seq1: str, seq2: str, temp: float = DEFAULT_TEMP) -> float:

    RNA.cvar.temperature = temp

    return RNA.duplexfold(

        normalize_rna(seq1),

        normalize_rna(seq2),

    ).energy





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





def window_scan(

    seq: str,

    abs_offset: int,

    k: int,

    temp: float,

) -> list[WindowHit]:

    seq = normalize_rna(seq)



    if k <= 0:

        raise ValueError("window_k must be > 0.")

    if len(seq) < k:

        return []



    hits: list[WindowHit] = []



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





def merge_top_hits(

    hits: list[WindowHit],

    source_seq: str,

    source_start: int,

    temp: float,

    tolerance: float = DG_TOL,

) -> list[MergedHit]:

    if not hits:

        return []



    best_dg = hits[0].dg

    top_hits = [h for h in hits if abs(h.dg - best_dg) <= tolerance]

    top_hits.sort(key=lambda h: h.abs_start)



    merged_intervals: list[tuple[int, int]] = []

    cur_start = top_hits[0].abs_start

    cur_end = top_hits[0].abs_end



    for h in top_hits[1:]:

        if h.abs_start <= cur_end + 1:

            cur_end = max(cur_end, h.abs_end)

        else:

            merged_intervals.append((cur_start, cur_end))

            cur_start, cur_end = h.abs_start, h.abs_end



    merged_intervals.append((cur_start, cur_end))



    merged: list[MergedHit] = []



    for start, end in merged_intervals:

        rel_start = start - source_start

        rel_end = end - source_start + 1

        subseq = source_seq[rel_start:rel_end]



        merged.append(

            MergedHit(

                abs_start=start,

                abs_end=end,

                subseq=subseq,

                dg=duplex_dg(subseq, revcomp_rna(subseq), temp),

            )

        )



    merged.sort(key=lambda x: x.dg)

    return merged





def build_rbs_from_scan_window(

    scan_seq: str,

    scan_abs_start: int,

    core_hit: WindowHit,

) -> tuple[str, str, str, str]:

    """

    Build the full o-RBS seed from the entire scanned ASD window.



    The strongest WindowHit defines the protected core, but the complete

    reverse-complemented scan window defines the full RBS. This guarantees that

    sequence surrounding the core is retained and can be mutated downstream.



    Returns:

        full_rbs, rbs_left, rbs_core, rbs_right

    """

    scan_seq = normalize_rna(scan_seq)

    full_rbs = revcomp_rna(scan_seq)



    # Coordinates of the core inside the original scan window.

    rel_core_start = core_hit.abs_start - scan_abs_start

    rel_core_end_exclusive = core_hit.abs_end - scan_abs_start + 1



    if not (

        0 <= rel_core_start < rel_core_end_exclusive <= len(scan_seq)

    ):

        raise RuntimeError(

            "ORBS core coordinates fall outside the selected ASD scan window."

        )



    # Reverse-complementing flips the coordinates.

    core_start_in_rbs = len(scan_seq) - rel_core_end_exclusive

    core_end_in_rbs = len(scan_seq) - rel_core_start



    rbs_left = full_rbs[:core_start_in_rbs]

    rbs_core = full_rbs[core_start_in_rbs:core_end_in_rbs]

    rbs_right = full_rbs[core_end_in_rbs:]



    expected_core = revcomp_rna(core_hit.subseq)

    if rbs_core != expected_core:

        raise RuntimeError(

            "Internal ORBS coordinate error: extracted RBS core does not match "

            "the reverse complement of the best hit."

        )



    return full_rbs, rbs_left, rbs_core, rbs_right





def build_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(

        description=(

            "Identify an orthogonal RBS core and export a full RBS seed "

            "containing mutable sequence around that protected core."

        )

    )



    parser.add_argument(

        "--asd-start",

        type=int,

        help="1-based start of the ASD search region.",

    )

    parser.add_argument(

        "--asd-end",

        type=int,

        help="1-based end of the ASD search region.",

    )

    parser.add_argument(

        "--output-dir",

        default=".",

        help="Directory where the <gene>_rbs_core.json file is written.",

    )



    # o-RBS / anti-SD input.

    orbs_source = parser.add_mutually_exclusive_group(required=True)

    orbs_source.add_argument(

        "--fasta",

        help="Anti-SD / o-RBS search sequence as FASTA.",

    )

    orbs_source.add_argument(

        "--sequence",

        help="Anti-SD / o-RBS search sequence as raw RNA/DNA.",

    )



    # Fixed 5' UTR.

    mrna5_source = parser.add_mutually_exclusive_group(required=True)

    mrna5_source.add_argument(

        "--mrna5",

        help="Fixed 5' UTR upstream of the standby/RBS; raw sequence.",

    )

    mrna5_source.add_argument(

        "--mrna5-fasta",

        help="Fixed 5' UTR upstream of the standby/RBS; FASTA file.",

    )



    parser.add_argument(

        "--standby",

        default="",

        help="Mutable standby sequence between the fixed 5' UTR and RBS.",

    )



    # Fixed CDS start.

    mrna3_source = parser.add_mutually_exclusive_group(required=True)

    mrna3_source.add_argument(

        "--mrna3",

        help="CDS start downstream of the spacer; raw sequence.",

    )

    mrna3_source.add_argument(

        "--mrna3-fasta",

        help="CDS start downstream of the spacer; FASTA file.",

    )



    # Retained for compatibility / reference in the output.

    crbs_source = parser.add_mutually_exclusive_group()

    crbs_source.add_argument(

        "--canonical-rbs",

        default="AUUCCUCCACUAG",

        help="Canonical RBS reference sequence; raw sequence.",

    )

    crbs_source.add_argument(

        "--canonical-rbs-fasta",

        help="Canonical RBS reference sequence; FASTA file.",

    )



    parser.add_argument(

        "--asd-tail-nt",

        type=int,

        default=DEFAULT_ASD_TAIL_NT,

        help=(

            "Number of nucleotides from the 3' end of the anti-SD input to use "

            "as the full RBS seed window."

        ),

    )

    parser.add_argument(

        "--window-k",

        type=int,

        default=DEFAULT_WINDOW_K,

        help="Window length used to identify the protected RBS core.",

    )

    parser.add_argument(

        "--temp",

        type=float,

        default=DEFAULT_TEMP,

    )

    parser.add_argument(

        "--name",

        default=DEFAULT_NAME,

        help="Gene/design name used for output naming.",

    )



    return parser





def main() -> None:

    args = build_parser().parse_args()



    # Resolve anti-SD / ORBS search sequence.

    search_sequence, inferred_name = resolve_input_sequence(

        args.sequence,

        args.fasta,

    )



    # Resolve fixed 5' UTR.

    if args.mrna5_fasta:

        _, five_prime_utr = load_fasta_sequence(

            args.mrna5_fasta,

            to_rna=True,

        )

    else:

        five_prime_utr = normalize_rna(args.mrna5)



    # Resolve fixed CDS start.

    if args.mrna3_fasta:

        _, cds_start = load_fasta_sequence(

            args.mrna3_fasta,

            to_rna=True,

        )

    else:

        cds_start = normalize_rna(args.mrna3)



    # Resolve canonical RBS reference.

    if args.canonical_rbs_fasta:

        _, canonical_rbs = load_fasta_sequence(

            args.canonical_rbs_fasta,

            to_rna=True,

        )

    else:

        canonical_rbs = normalize_rna(args.canonical_rbs)



    standby = normalize_rna(args.standby)

    gene_name = sanitize_name(args.name or inferred_name)



    if not search_sequence:

        raise SystemExit("The ORBS/anti-SD search sequence is empty.")



    # Define the ASD scan window. Clamp the default tail length to the actual

    # sequence length so a short anti-SD sequence works correctly.

    end = args.asd_end or len(search_sequence)



    if args.asd_start is not None:

        start = args.asd_start

    else:

        tail_nt = max(1, min(int(args.asd_tail_nt), len(search_sequence)))

        start = len(search_sequence) - tail_nt + 1



    if start < 1 or end > len(search_sequence) or start > end:

        raise SystemExit("Invalid ASD window coordinates.")



    scan_seq = search_sequence[start - 1:end]

    scan_abs_start = start



    if len(scan_seq) < int(args.window_k):

        raise SystemExit(

            f"ASD scan window is {len(scan_seq)} nt, shorter than "

            f"window_k={args.window_k}."

        )



    hits = window_scan(

        scan_seq,

        scan_abs_start,

        int(args.window_k),

        float(args.temp),

    )



    print("RAW WINDOW HITS")

    print("ABS_START\tABS_END\tSUBSEQ\tDG")

    for h in hits:

        print(f"{h.abs_start}\t{h.abs_end}\t{h.subseq}\t{h.dg:.2f}")



    merged = merge_top_hits(

        hits,

        scan_seq,

        scan_abs_start,

        float(args.temp),

    )



    print("MERGED TOP-HIT REGIONS")

    print("ABS_START\tABS_END\tLEN\tSEQUENCE\tDG")

    for m in merged:

        print(

            f"{m.abs_start}\t{m.abs_end}\t{len(m.subseq)}\t"

            f"{m.subseq}\t{m.dg:.2f}"

        )



    if not hits:

        print("No ORBS core hit found; no JSON file written.")

        return



    # IMPORTANT CHANGE:

    # The best hit defines the protected core, while the entire ASD scan window

    # becomes the full RBS seed. Previously full_rbs came from merged[0].subseq,

    # which could be identical to the core and leave rbs_left/rbs_right empty.

    core_hit = hits[0]

    full_rbs, rbs_left, rbs_core, rbs_right = build_rbs_from_scan_window(

        scan_seq,

        scan_abs_start,

        core_hit,

    )



    print("SELECTED FULL RBS SEED")

    print(f"ASD scan window : {scan_seq}")

    print(f"Full RBS        : {full_rbs}")

    print(f"RBS left        : {rbs_left or '<empty>'}")

    print(f"RBS core        : {rbs_core}")

    print(f"RBS right       : {rbs_right or '<empty>'}")



    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{gene_name}_rbs_core.json"



    candidate = {

        "name": gene_name,

        "five_prime_utr": five_prime_utr,

        "standby": standby,

        "canonical_rbs": canonical_rbs,

        "rbs_left": rbs_left,

        "rbs_core": rbs_core,

        "rbs_right": rbs_right,

        "rbs": full_rbs,

        "spacer": "",

        "cds_start": cds_start,

        "mutable_regions": [

            "standby",

            "rbs_left",

            "rbs_right",

            "spacer",

        ],

        "protected_regions": [

            "five_prime_utr",

            "rbs_core",

            "cds_start",

        ],

    }



    with output_path.open("w") as fh:

        json.dump([candidate], fh, indent=4)

        fh.write("\n")



    print(f"WROTE {output_path}")





if __name__ == "__main__":

    main()
