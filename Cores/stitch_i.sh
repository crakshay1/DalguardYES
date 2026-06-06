#!/bin/sh

python3 orbs_duplex.py --fasta Tv3test.fa

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"

RNA5_SEQ="tacacacgaataaaagataacaaagatgagtaaaggagaagaactttt"
RNA3_SEQ="atgagtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggcgatgttaatgggcaaaaattctctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactgggaagctacctgttccatggccaacacttgtcactactttctcttatggtgttcaatgcttttcaagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttttacaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaaatggaatacaactataactcacataatgtatatacatcatggcagacaaaccaaagaatggaatcaaagttaacttcaaaattagacacaacattaaagatggaagcgttcaattagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttccaaagatcccaacgaaaagagagatcacatgatccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacaaataaatgtccagacttccaattgacactaaagtgtccgaacaattactaaattctcagggttcctggttaaattcaggctgagactttatttatatatttatagattcattaaaattttatgaataatttattgatgttattataaggggctatttttcttattaaatagggctactggagtgtat"
ORBS="oRibo_rbs_core.fa"

RNA5_FILE=$(mktemp "${TMPDIR:-/tmp}/mrna5.XXXXXX.fa")
RNA3_FILE=$(mktemp "${TMPDIR:-/tmp}/mrna3.XXXXXX.fa")
trap 'rm -f "$RNA5_FILE" "$RNA3_FILE"' EXIT

cat > "$RNA5_FILE" <<EOF
>mrna5
$RNA5_SEQ
EOF

cat > "$RNA3_FILE" <<EOF
>mrna3
$RNA3_SEQ
EOF

python3 mrna_stitch.py --mrna5 "$RNA5_FILE" --mrna3 "$RNA3_FILE" --orbs "$ORBS" --name gfp
