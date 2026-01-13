#include <stdio.h>
#include <string.h>
#include "heat.h"
#include "grid.h"
#include "config.h"
// SIMDe SSE (single-precision) for portable `simde__m128` operations
#include <simde/x86/sse.h>

int heat_step_seq(grid_t *g, grid_t *g_copy)
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

    // actualizar solo los puntos interiores (vectorizado con SIMDe `simde__m128`)
    for (int y = 1; y < n_filas - 1; ++y)
    {
        int x = 1;
        const int vec_width = 4; // 4 floats per `simde__m128`
        // For lanes x..x+3, right neighbor accesses up to x+4, which must be <= n_cols-1
        int x_end = n_cols - 5;

        if (x_end >= x)
        {
            simde__m128 coeff = simde_mm_set1_ps((real_t)0.25f);
            for (; x <= x_end; x += vec_width)
            {
                const real_t *up_ptr = &GRID_AT(g, y - 1, x);
                const real_t *down_ptr = &GRID_AT(g, y + 1, x);
                const real_t *left_ptr = &GRID_AT(g, y, x - 1);
                const real_t *right_ptr = &GRID_AT(g, y, x + 1);

                simde__m128 v_up = simde_mm_loadu_ps(up_ptr);
                simde__m128 v_down = simde_mm_loadu_ps(down_ptr);
                simde__m128 v_left = simde_mm_loadu_ps(left_ptr);
                simde__m128 v_right = simde_mm_loadu_ps(right_ptr);

                simde__m128 sum = simde_mm_add_ps(v_up, v_down);
                sum = simde_mm_add_ps(sum, v_left);
                sum = simde_mm_add_ps(sum, v_right);

                simde__m128 res = simde_mm_mul_ps(sum, coeff);

                real_t *dst = &GRID_AT(g_copy, y, x);
                simde_mm_storeu_ps(dst, res);
            }
        }

        // scalar fallback for remaining columns
        for (; x < n_cols - 1; ++x)
        {
            real_t up = GRID_AT(g, y - 1, x);
            real_t down = GRID_AT(g, y + 1, x);
            real_t left = GRID_AT(g, y, x - 1);
            real_t right = GRID_AT(g, y, x + 1);

            GRID_AT(g_copy, y, x) = (up + down + left + right) * (real_t)0.25f;
        }
    }

    // copiar malla temporal a malla original
    grid_swap(g, g_copy);
    return 0;
}

int heat_iterate_seq(grid_t *g, int nsteps)
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
        check_res = heat_step_seq(g, &g_copy);
    }
    grid_free(&g_copy);
    return check_res;
}
