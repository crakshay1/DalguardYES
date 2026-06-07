#!/usr/bin/env python3
"""
Usage:
    python top_candidate.py response.json
    curl ... | python top_candidate.py
"""

import json
import sys


def to_fasta(header: str, seq: str, width: int = 80) -> str:
    lines = [f">{header}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i:i + width])
    return "\n".join(lines)


def main() -> None:
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            data = json.load(f)
    else:
        data = json.load(sys.stdin)

    candidates = data.get("candidates", [])
    if not candidates:
        sys.exit("No candidates in response.")

    top = candidates[0]
    seq = top["fullSeq"].upper().replace("T", "U")

    header = (
        f"{top['candidateId']} "
        f"fitness={top['fitness']} "
        f"status={top.get('status', 'unknown')} "
        f"orthTIR={float(top['orthTIR']):.4g} "
        f"wtTIR={float(top['wtTIR']):.4g} "
        f"tScore={float(top['tScore']):.4g}"
    )

    print(to_fasta(header, seq))


if __name__ == "__main__":
    main()
