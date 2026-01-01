#include <stdio.h>
#include "grid.h"
#include "config.h"

int heat_step(grid_t *g)
{
    if (!g || !g->data)
        return -1;

    int n_cols = g->nx;
    int n_filas = g->ny;

    // reservar malla temporal
    grid_t *g_copy;
    if (grid_alloc(g_copy, n_cols, n_filas) != 0)
    {
        return -2;
    }

    int max = n_cols * n_filas;
    real_t *data = g->data;
    for (int count = 0; count < max; count++)
    {
        int x = count % n_cols;
        int y = count / n_cols;

        int n, sum;
        int vals[4];

        real_t mean;
        ;
        // operar para bordes
        if (y > 0 && y < n_cols - 1 && x > 0 && x < n_filas - 1)
        {
            n = 0;
            vals[n++] = GRID_AT(g, y - 1, x); // arriba
            vals[n++] = GRID_AT(g, y + 1, x); // abajo
            vals[n++] = GRID_AT(g, y, x - 1); // izquierda
            vals[n++] = GRID_AT(g, y, x + 1); // derecha

        }
        else if (y == 0) // si y == 0 -> borde superior
        {
            n = 0;
            vals[n++] = GRID_AT(g, y + 1, x); // abajo
            vals[n++] = GRID_AT(g, y, x - 1); // izquierda
            vals[n++] = GRID_AT(g, y, x + 1); // derecha
        }
        else if (y == n_cols - 1) // si y == n_cols - 1 -> borde inferior
        {
            n = 0;
            vals[n++] = GRID_AT(g, y - 1, x); // arriba
            vals[n++] = GRID_AT(g, y, x - 1); // izquierda
            vals[n++] = GRID_AT(g, y, x + 1); // derecha
        }
        else if (x == 0) // si x == 0 -> borde izquierdo
        {
            n = 0;
            vals[n++] = GRID_AT(g, y - 1, x); // arriba
            vals[n++] = GRID_AT(g, y + 1, x); // abajo
            vals[n++] = GRID_AT(g, y, x + 1); // derecha
        }
        else if (x == n_filas - 1) // si x == n_filas - 1 -> borde derecho
        {
            n = 0;
            vals[n++] = GRID_AT(g, y - 1, x); // arriba
            vals[n++] = GRID_AT(g, y + 1, x); // abajo
            vals[n++] = GRID_AT(g, y, x - 1); // izquierda
        }

        for (size_t i = 0; i < n; i++)
        {
            sum += vals[i];
        }
        mean = (real_t) sum / n; //cast to real_t so it gives a real_t result
        GRID_AT(g_copy, y, x) = mean;
    }

    // copiar malla temporal a malla original
    grid_swap(g, g_copy);
    grid_free(g_copy);
}

int heat_iterate(grid_t *g, int nsteps)
{
    int check_res;
    for (size_t i = 0; i < nsteps && check_res == 0; i++)
    {
        check_res = heat_step(g);
    }
}
