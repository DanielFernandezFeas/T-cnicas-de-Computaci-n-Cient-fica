#include "grid.h"
#include "config.h"
#include <stdio.h>
#include <stdlib.h>

int grid_alloc(grid_t *g, int nx, int ny) {
    if (!g || nx <= 0 || ny <= 0) return -1;

    g->nx = nx;
    g->ny = ny;

    size_t n = (size_t)nx * (size_t)ny;
    g->data = (real_t *)malloc(n * sizeof(real_t));
    if (!g->data) return -2;

    return 0;
}

void grid_free(grid_t *g) {
    if (!g) return;
    free(g->data);
    g->data = NULL;
    g->nx = 0;
    g->ny = 0;
}

void grid_fill(grid_t *g, real_t v) {
    if (!g || !g->data) return;

    size_t n = (size_t)g->nx * (size_t)g->ny;
    for (size_t i = 0; i < n; ++i) {
        g->data[i] = v;
    }
}

/*
 * Aplica condiciones de contorno (Boundary Conditions, BC).
 *
 * Convención que estamos usando (según lo visto en tu salida):
 * - Borde superior (y=0): hot
 * - Bordes izquierdo (x=0), derecho (x=nx-1) e inferior (y=ny-1): cold
 *
 * Nota: Fíjate que esto hace que las ESQUINAS superiores queden a cold,
 * porque después de poner el borde superior a hot, forzamos x=0 y x=nx-1 a cold.
 * Por eso en tu salida el checksum sale 800 en vez de 1000 para una malla 10x10.
 */
void grid_apply_bc(grid_t *g, real_t hot, real_t cold) {
    if (!g || !g->data) return;

    // Borde superior a hot
    for (int x = 0; x < g->nx; ++x) {
        GRID_AT(g, 0, x) = hot;
    }

    // Bordes restantes a cold (izq, der, inferior)
    for (int y = 0; y < g->ny; ++y) {
        GRID_AT(g, y, 0) = cold;
        GRID_AT(g, y, g->nx - 1) = cold;
    }

    for (int x = 0; x < g->nx; ++x) {
        GRID_AT(g, g->ny - 1, x) = cold;
    }
}

void grid_init_bc(grid_t *g, real_t hot, real_t cold) {
    // Primero todo cold, luego aplicamos BC
    grid_fill(g, cold);
    grid_apply_bc(g, hot, cold);
}

void grid_swap(grid_t *a, grid_t *b) {
    if (!a || !b) return;

    if (a->nx != b->nx || a->ny != b->ny) return;

    real_t *tmp = a->data;
    a->data = b->data;
    b->data = tmp;
}

real_t grid_checksum(const grid_t *g) {
    if (!g || !g->data) return (real_t)0;

    size_t n = (size_t)g->nx * (size_t)g->ny;
    real_t sum = 0;
    for (size_t i = 0; i < n; ++i) {
        sum += g->data[i];
    }
    return sum;
}

void grid_print(const grid_t *g, int max_nx, int max_ny) {
    if (!g || !g->data) return;

    int ny = g->ny < max_ny ? g->ny : max_ny;
    int nx = g->nx < max_nx ? g->nx : max_nx;

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            printf("%7.2f ", (double)GRID_AT(g, y, x));
        }
        if (nx < g->nx) printf(" ...");
        printf("\n");
    }
    if (ny < g->ny) printf("...\n");
}
