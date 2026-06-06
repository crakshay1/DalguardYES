#!/usr/bin/env python3
"""
mrna_stitch.py

Stitch an mRNA fragment upstream of the RBS, the o-RBS sequence, and the
downstream mRNA fragment into one recombinant transcript and write it as FASTA.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re


def normalize_rna(seq: str) -> str:
	return seq.upper().replace("T", "U")


def sanitize_name(name: str) -> str:
	cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
	return cleaned.strip("._-") or "sequence"


def load_fasta_sequence(path: str) -> tuple[str, str]:
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
				parts.append(normalize_rna(line))
	return header, "".join(parts)


def load_first_fasta_record_sequence(path: str) -> tuple[str, str]:
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
				continue
			if header:
				parts.append(normalize_rna(line))
	return header, "".join(parts)


def write_fasta(path: Path, header: str, sequence: str, width: int = 80) -> None:
	sequence = normalize_rna(sequence)
	with path.open("w") as fh:
		fh.write(f">{header}\n")
		for i in range(0, len(sequence), width):
			fh.write(sequence[i:i + width] + "\n")


def resolve_sequence(value: str) -> tuple[str, str | None]:
	candidate = Path(value)
	if candidate.is_file():
		header, sequence = load_first_fasta_record_sequence(value)
		return sequence, header
	return normalize_rna(value), None


def infer_gene_name(explicit_name: str | None, orbs_header: str | None, mrna5_header: str | None, mrna3_header: str | None) -> str:
	if explicit_name:
		return sanitize_name(explicit_name)

	for header in (orbs_header, mrna5_header, mrna3_header):
		if not header:
			continue
		token = header.split()[0]
		token = token[2:] if token.startswith("o_") else token
		token = token[:-3] if token.endswith("_cr") else token
		if token:
			return sanitize_name(token)

	raise SystemExit("--name is required when all inputs are manual sequences.")


def build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(description="Stitch mRNA fragments with an o-RBS and export the recombinant transcript as FASTA.")
	parser.add_argument("--mrna5", required=True, help="mRNA fragment upstream of the RBS; FASTA path or raw sequence.")
	parser.add_argument("--orbs", required=True, help="o-RBS sequence produced by asd_duplex.py; FASTA path or raw sequence.")
	parser.add_argument("--mrna3", required=True, help="mRNA fragment downstream of the RBS; FASTA path or raw sequence.")
	parser.add_argument("--name", help="Gene name used for the output FASTA file.")
	parser.add_argument("--output-dir", default=".", help="Directory where the stitched FASTA file is written.")
	return parser


def main() -> None:
	args = build_parser().parse_args()

	mrna5_seq, mrna5_header = resolve_sequence(args.mrna5)
	orbs_seq, orbs_header = resolve_sequence(args.orbs)
	mrna3_seq, mrna3_header = resolve_sequence(args.mrna3)

	gene_name = infer_gene_name(args.name, orbs_header, mrna5_header, mrna3_header)
	stitched = mrna5_seq + orbs_seq + mrna3_seq

	output_dir = Path(args.output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	output_path = output_dir / f"{gene_name}_cr.fa"

	write_fasta(output_path, f"{gene_name}_cr", stitched)
	print(f"WROTE {output_path}")


if __name__ == "__main__":
	main()
