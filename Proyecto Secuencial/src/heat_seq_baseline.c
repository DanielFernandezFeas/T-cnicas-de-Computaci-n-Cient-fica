#include <string.h>
#include "heat.h"
#include "grid.h"
#include "config.h"

/*
 * BASELINE (escalares):
 * - Copia toda la malla con memcpy para preservar bordes/BC
 * - Actualiza solo el interior con bucles escalares
 * - Hace swap al final de cada paso
 */

int heat_step_seq(grid_t *g, grid_t *g_copy)
{
    if (!g || !g->data) return -1;
    if (!g_copy) return -2;

    int nx = g->nx;
    int ny = g->ny;

    // Asegurar que g_copy tiene memoria y tamaño correcto
    if (g_copy->data == NULL || g_copy->nx != nx || g_copy->ny != ny) {
        // Si tenía memoria previa pero tamaño distinto, libérala antes
        if (g_copy->data != NULL) {
            grid_free(g_copy);
        }
        if (grid_alloc(g_copy, nx, ny) != 0) return -2;
    }

    // Copia completa: mantiene bordes/condiciones de contorno
    size_t n = (size_t)nx * (size_t)ny;
    memcpy(g_copy->data, g->data, n * sizeof(real_t));

    // Actualiza interior (Jacobi 2D)
    for (int y = 1; y < ny - 1; ++y) {
        for (int x = 1; x < nx - 1; ++x) {
            real_t up    = GRID_AT(g, y - 1, x);
            real_t down  = GRID_AT(g, y + 1, x);
            real_t left  = GRID_AT(g, y, x - 1);
            real_t right = GRID_AT(g, y, x + 1);

            GRID_AT(g_copy, y, x) = (up + down + left + right) * (real_t)0.25;
        }
    }

    // Swap: g pasa a ser el nuevo estado
    grid_swap(g, g_copy);
    return 0;
}

int heat_iterate_seq(grid_t *g, int nsteps)
{
    if (!g || !g->data) return -1;
    if (nsteps <= 0) return 0;

    int rc = 0;
    grid_t tmp = (grid_t){0};

    // Reserva una vez el grid temporal
    if (grid_alloc(&tmp, g->nx, g->ny) != 0) return -2;

    for (int i = 0; i < nsteps && rc == 0; ++i) {
        rc = heat_step_seq(g, &tmp);
    }

    grid_free(&tmp);
    return rc;
}
