#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/raw

CSV="results/raw/omp.csv"
ITERS=1000
REPS=5

echo "version,nx,ny,iters,threads_or_procs,time_s,checksum" > "$CSV"

for N in 100 500 1000; do
  for T in 1 2 4 8 16; do
    ./main_omp "$N" "$N" "$ITERS" "$REPS" "$T" >> "$CSV"
  done
done

echo "OK -> generado $CSV"
