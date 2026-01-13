#include "heat.h"
#include "grid.h"
#include "config.h"

/* Copia SOLO los bordes para mantener BC */
static void copy_boundaries(const grid_t *src, grid_t *dst) {
    int nx = src->nx, ny = src->ny;

    for (int x = 0; x < nx; ++x) {
        GRID_AT(dst, 0, x)      = GRID_AT(src, 0, x);
        GRID_AT(dst, ny - 1, x) = GRID_AT(src, ny - 1, x);
    }
    for (int y = 0; y < ny; ++y) {
        GRID_AT(dst, y, 0)      = GRID_AT(src, y, 0);
        GRID_AT(dst, y, nx - 1) = GRID_AT(src, y, nx - 1);
    }
}

/* MISMO nombre que tu main usa: heat_iterate_seq */
int heat_iterate_seq(grid_t *g, int nsteps) {
    if (!g || !g->data) return -1;
    if (nsteps <= 0) return 0;

    grid_t tmp = (grid_t){0};
    if (grid_alloc(&tmp, g->nx, g->ny) != 0) return -2;

    int nx = g->nx, ny = g->ny;

    for (int i = 0; i < nsteps; ++i) {
        copy_boundaries(g, &tmp);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                real_t up    = GRID_AT(g, y - 1, x);
                real_t down  = GRID_AT(g, y + 1, x);
                real_t left  = GRID_AT(g, y, x - 1);
                real_t right = GRID_AT(g, y, x + 1);
                GRID_AT(&tmp, y, x) = (up + down + left + right) * (real_t)0.25;
            }
        }

        grid_swap(g, &tmp);
    }

    grid_free(&tmp);
    return 0;
}
