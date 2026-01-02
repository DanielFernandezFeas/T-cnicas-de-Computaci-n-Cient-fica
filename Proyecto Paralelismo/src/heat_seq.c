#include <stdio.h>
#include <string.h>
#include "heat.h"
#include "grid.h"
#include "config.h"

int heat_step(grid_t *g, grid_t *g_copy)
{
    if (!g || !g->data)
        return -1;

    int n_cols = g->nx;
    int n_filas = g->ny;

    // require a valid temporary grid pointer
    if (g_copy == NULL)
    {
        return -2;
    }

    // allocate temp grid if its storage is missing or size mismatches
    if (g_copy->data == NULL || g_copy->nx != n_cols || g_copy->ny != n_filas)
    {
        if (grid_alloc(g_copy, n_cols, n_filas) != 0)
        {
            return -2;
        }
    }


    size_t max = (size_t)n_cols * (size_t)n_filas;
    memcpy(g_copy->data, g->data, max * sizeof(real_t));

    // actualizar solo los puntos interiores
    for (int y = 1; y < n_filas - 1; ++y)
    {
        for (int x = 1; x < n_cols - 1; ++x)
        {
            real_t up = GRID_AT(g, y - 1, x);
            real_t down = GRID_AT(g, y + 1, x);
            real_t left = GRID_AT(g, y, x - 1);
            real_t right = GRID_AT(g, y, x + 1);

            GRID_AT(g_copy, y, x) = (up + down + left + right) * (real_t)0.25;
        }
    }

    // copiar malla temporal a malla original
    grid_swap(g, g_copy);
    return 0;
}

int heat_iterate(grid_t *g, int nsteps)
{
    int check_res = 0;
    grid_t g_copy;
    int n_cols = g->nx;
    int n_filas = g->ny;
    if (grid_alloc(&g_copy, n_cols, n_filas) != 0)
    {
        return -2;
    }

    for (size_t i = 0; i < nsteps && check_res == 0; i++)
    {
        check_res = heat_step(g, &g_copy);
    }
    grid_free(&g_copy);
    return check_res;
}
