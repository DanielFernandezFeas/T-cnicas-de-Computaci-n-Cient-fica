/*
 * Simplified test harness: accepts three command-line arguments:
 *   ./main_seq <nx> <ny> <iters> [reps]
 * where reps is optional (default 1). Produces a CSV line:
 * seq,nx,ny,iters,procs,time_s,checksum
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "grid.h"
#include "config.h"
#include "heat.h"

static double now_sec(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <nx> <ny> <iters> [reps]\n", argv[0]);
        return 1;
    }

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int iters = atoi(argv[3]);
    int reps = (argc > 4) ? atoi(argv[4]) : 1;

    if (nx < 3 || ny < 3 || iters <= 0 || reps <= 0) {
        fprintf(stderr, "Invalid args: nx,ny >= 3; iters,reps > 0\n");
        return 1;
    }

    const real_t hot = (real_t)HOT;
    const real_t cold = (real_t)COLD;

    grid_t g;
    if (grid_alloc(&g, nx, ny) != 0) {
        fprintf(stderr, "grid_alloc failed\n");
        return 1;
    }

    double best = 1e300;
    double sum = 0.0;
    real_t checksum_last = 0;

    for (int r = 0; r < reps; ++r) {
        grid_init_bc(&g, hot, cold);

        double t0 = now_sec();
        int rc = heat_iterate_seq(&g, iters);
        double t1 = now_sec();

        if (rc != 0) {
            fprintf(stderr, "heat_iterate_seq failed (rc=%d)\n", rc);
            grid_free(&g);
            return 1;
        }

        double dt = t1 - t0;
        sum += dt;
        if (dt < best) best = dt;

        checksum_last = grid_checksum(&g);
    }

    double avg = sum / reps;

    // CSV output for automated collection
    printf("seq,%d,%d,%d,%d,%.9f,%.6f\n", nx, ny, iters, 1, avg, (double)checksum_last);
    fprintf(stderr, "SEQ nx=%d ny=%d iters=%d reps=%d avg=%.6fs best=%.6fs checksum=%.6f\n",
            nx, ny, iters, reps, avg, best, (double)checksum_last);

    grid_free(&g);
    return 0;
}
