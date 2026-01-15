#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "grid.h"
#include "config.h"
#include "heat.h"

static double now_sec(void) {
    return omp_get_wtime();
}


// Uso: ./main_omp nx ny iters reps threads
int main(int argc, char **argv) {
    int nx      = (argc > 1) ? atoi(argv[1]) : 100;
    int ny      = (argc > 2) ? atoi(argv[2]) : 100;
    int iters   = (argc > 3) ? atoi(argv[3]) : 1000;
    int reps    = (argc > 4) ? atoi(argv[4]) : 5;
    int threads = (argc > 5) ? atoi(argv[5]) : 1;

    if (nx < 3 || ny < 3 || iters <= 0 || reps <= 0 || threads <= 0) {
        fprintf(stderr, "Uso: %s nx ny iters reps threads\n", argv[0]);
        return 1;
    }

    const real_t hot  = (real_t)100.0;
    const real_t cold = (real_t)0.0;

    grid_t g;
    if (grid_alloc(&g, nx, ny) != 0) {
        fprintf(stderr, "Error reservando malla\n");
        return 1;
    }

    double sum = 0.0;
    double best = 1e300;
    real_t checksum_last = 0;

    for (int r = 0; r < reps; ++r) {
        grid_init_bc(&g, hot, cold);

        double t0 = now_sec();
        int rc = heat_iterate_omp(&g, iters, threads);
        double t1 = now_sec();

        if (rc != 0) {
            fprintf(stderr, "Error en heat_iterate_omp (rc=%d)\n", rc);
            grid_free(&g);
            return 1;
        }

        double dt = t1 - t0;
        sum += dt;
        if (dt < best) best = dt;

        checksum_last = grid_checksum(&g);
    }

    double avg = sum / reps;

    // CSV: version,nx,ny,iters,threads_or_procs,time_s,checksum
    printf("omp,%d,%d,%d,%d,%.9f,%.6f\n", nx, ny, iters, threads, avg, (double)checksum_last);

    fprintf(stderr, "OMP nx=%d ny=%d iters=%d reps=%d threads=%d avg=%.6fs best=%.6fs checksum=%.2f\n",
            nx, ny, iters, reps, threads, avg, best, (double)checksum_last);

    grid_free(&g);
    return 0;
}
