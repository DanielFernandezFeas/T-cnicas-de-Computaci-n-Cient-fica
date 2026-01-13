//
// Created by DANIEL on 30/12/2025.
//

#ifndef T_CNICAS_DE_COMPUTACI_N_CIENT_FICA_GRID_H
#define T_CNICAS_DE_COMPUTACI_N_CIENT_FICA_GRID_H

#include <stddef.h>
#include "config.h"

typedef struct {
    int nx, ny;
    real_t *data;
} grid_t;

#define GRID_IDX(g, y, x) ((size_t)(y) * (size_t)((g)->nx) + (size_t)(x))
#define GRID_AT(g, y, x)  ((g)->data[GRID_IDX((g),(y),(x))])

int  grid_alloc(grid_t *g, int nx, int ny);
void grid_free(grid_t *g);

void grid_fill(grid_t *g, real_t v);
void grid_apply_bc(grid_t *g, real_t hot, real_t cold);
void grid_init_bc(grid_t *g, real_t hot, real_t cold);
void grid_print(const grid_t *g, int max_nx, int max_ny);


void grid_swap(grid_t *a, grid_t *b);
real_t grid_checksum(const grid_t *g);

#endif //T_CNICAS_DE_COMPUTACI_N_CIENT_FICA_GRID_H