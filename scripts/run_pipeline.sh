#!/usr/bin/env bash
set -euo pipefail

INPUT="data/somatic_mutations.tsv"
OUTPUT="results/candidate_neoantigens.tsv"

echo "=== Neoantigen Extraction & Profiling Pipeline ==="
echo "[1/2] Processing somatic variants and extracting 9-mers..."
python scripts/extract_neoantigens.py "$INPUT" "$OUTPUT"

echo "[2/2] Top Candidate Peptides:"
echo "------------------------------------------------------------"
column -t -s $'\t' "$OUTPUT"
echo "------------------------------------------------------------"
echo "=== Pipeline Completed Successfully ==="
