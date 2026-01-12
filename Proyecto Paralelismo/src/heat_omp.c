#include <omp.h>
#include "heat.h"
#include "grid.h"
#include "config.h"

/* Copia SOLO los bordes de src -> dst para mantener las BC sin memcpy masivo */
static void copy_boundaries(const grid_t *src, grid_t *dst) {
    int nx = src->nx;
    int ny = src->ny;

    // fila superior e inferior
    for (int x = 0; x < nx; ++x) {
        GRID_AT(dst, 0, x)      = GRID_AT(src, 0, x);
        GRID_AT(dst, ny - 1, x) = GRID_AT(src, ny - 1, x);
    }

    // columna izquierda y derecha
    for (int y = 0; y < ny; ++y) {
        GRID_AT(dst, y, 0)      = GRID_AT(src, y, 0);
        GRID_AT(dst, y, nx - 1) = GRID_AT(src, y, nx - 1);
    }
}

/*
 * 1 iteración Jacobi en OpenMP:
 * - lee de g (old)
 * - escribe interior en g_copy (new)
 * - copia bordes (BC)
 * - swap(g, g_copy)
 */
int heat_step_omp(grid_t *g, grid_t *g_copy, int nthreads) {
    if (!g || !g->data) return -1;
    if (!g_copy) return -2;

    int nx = g->nx;
    int ny = g->ny;

    // asegurar que g_copy existe y tiene el mismo tamaño
    if (g_copy->data == NULL || g_copy->nx != nx || g_copy->ny != ny) {
        if (g_copy->data != NULL) {
            grid_free(g_copy);
        }
        if (grid_alloc(g_copy, nx, ny) != 0) return -2;
    }

    // mantener BC
    copy_boundaries(g, g_copy);

    // actualizar interior en paralelo (sin tocar bordes)
    if (nthreads > 0) {
        #pragma omp parallel for num_threads(nthreads) collapse(2) schedule(static)
        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                real_t up    = GRID_AT(g, y - 1, x);
                real_t down  = GRID_AT(g, y + 1, x);
                real_t left  = GRID_AT(g, y, x - 1);
                real_t right = GRID_AT(g, y, x + 1);

                GRID_AT(g_copy, y, x) = (up + down + left + right) * (real_t)0.25;
            }
        }
    } else {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                real_t up    = GRID_AT(g, y - 1, x);
                real_t down  = GRID_AT(g, y + 1, x);
                real_t left  = GRID_AT(g, y, x - 1);
                real_t right = GRID_AT(g, y, x + 1);

                GRID_AT(g_copy, y, x) = (up + down + left + right) * (real_t)0.25;
            }
        }
    }

    grid_swap(g, g_copy);
    return 0;
}

int heat_iterate_omp(grid_t *g, int nsteps, int nthreads) {
    if (!g || !g->data) return -1;
    if (nsteps <= 0) return 0;

    grid_t tmp = {0};
    if (grid_alloc(&tmp, g->nx, g->ny) != 0) return -2;

    int rc = 0;
    for (int i = 0; i < nsteps && rc == 0; ++i) {
        rc = heat_step_omp(g, &tmp, nthreads);
    }

    grid_free(&tmp);
    return rc;
}
