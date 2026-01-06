#!/usr/bin/env bash
# Usage: ./run_seq.sh [binary] [out_csv]
# Defaults: binary=../apps/main_seq  out_csv=../results/raw/seq.csv

set -euo pipefail

BINARY=${1:-../apps/main_seq}
OUT_CSV=${2:-../results/raw/seq.csv}

mkdir -p "$(dirname "$OUT_CSV")"

echo "# timestamp,csv_line" > "$OUT_CSV"

# sizes and iterations to test (modify as needed)
SIZES=(64 128 256 512)
ITERS=(100 500 1000)
REPS=1

for s in "${SIZES[@]}"; do
  for it in "${ITERS[@]}"; do
    echo "Running size=${s}x${s} iters=${it} reps=${REPS}" >&2
    # run the binary; it prints CSV line: seq,nx,ny,iters,procs,time_s,checksum
    CSV_LINE=$($BINARY "$s" "$s" "$it" "$REPS")
    printf "%s,%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$CSV_LINE" >> "$OUT_CSV"
  done
done

echo "Results written to $OUT_CSV" >&2
