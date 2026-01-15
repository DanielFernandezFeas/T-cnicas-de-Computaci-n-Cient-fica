#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#include "grid.h"
#include "config.h"
#include "heat.h"

static double now_sec(void) {
    return MPI_Wtime();
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Parámetros SOLO los lee rank 0 */
    int nx, ny, iters, reps;

    if (rank == 0) {
        nx    = (argc > 1) ? atoi(argv[1]) : 100;
        ny    = (argc > 2) ? atoi(argv[2]) : 100;
        iters = (argc > 3) ? atoi(argv[3]) : 1000;
        reps  = (argc > 4) ? atoi(argv[4]) : 5;

        if (nx < 3 || ny < 3 || iters <= 0 || reps <= 0) {
            fprintf(stderr, "Uso: %s nx ny iters reps\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Broadcast de parámetros */
    MPI_Bcast(&nx,    1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ny,    1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&reps,  1, MPI_INT, 0, MPI_COMM_WORLD);

    grid_t g = {0};

    /* SOLO rank 0 tiene la malla global */
    if (rank == 0) {
        if (grid_alloc(&g, nx, ny) != 0) {
            fprintf(stderr, "Error reservando malla\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Ensure dimensions are available on all ranks (data only on rank 0) */
    g.nx = nx;
    g.ny = ny;

    double sum = 0.0, best = 1e300;
    real_t checksum_last = 0;

    for (int r = 0; r < reps; ++r) {

        if (rank == 0) {
            grid_init_bc(&g, HOT, COLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        double t0 = now_sec();
        int rc = heat_iterate_mpi(&g, iters, MPI_COMM_WORLD);
        double t1 = now_sec();

        if (rc != 0) {
            if (rank == 0)
                fprintf(stderr, "Error en heat_iterate_mpi (rc=%d)\n", rc);
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        double dt_local = t1 - t0;
        double dt = 0.0;
        MPI_Reduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            sum += dt;
            if (dt < best) best = dt;
            checksum_last = grid_checksum(&g);
    }

    if (rank == 0) {
        double avg = sum / reps;

        printf("mpi,%d,%d,%d,%.9f,%.6f\n",
               nx, ny, iters, avg, (double)checksum_last);

        fprintf(stderr,
                "MPI nx=%d ny=%d iters=%d reps=%d procs=%d "
                "avg=%.6fs best=%.6fs checksum=%.2f\n",
                nx, ny, iters, reps, size,
                avg, best, (double)checksum_last);

        grid_free(&g);
    }

    MPI_Finalize();
    return 0;
}
