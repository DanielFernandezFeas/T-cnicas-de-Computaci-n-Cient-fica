#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/raw

CSV="results/raw/mpi.csv"
ERR_LOG="results/raw/mpi_errors.log"
ITERS=1000
REPS=5

# Cabecera CSV
echo "version,nx,ny,iters,threads_or_procs,time_s,checksum" > "$CSV"

# Limpiar log de errores previo
> "$ERR_LOG"

# Tamaños de malla
for N in 100 500 1000; do
    # Número de procesos MPI
    for P in 1 2 4 8 16; do
        echo "Ejecutando MPI: N=$N P=$P ..."
        
        # Ejecutar MPI
        mpirun -np "$P" ./main_mpi "$N" "$N" "$ITERS" "$REPS" \
            >> "$CSV" 2>> "$ERR_LOG"
    done
done

echo "OK -> generado $CSV"
echo "Errores (si los hay) en $ERR_LOG"
