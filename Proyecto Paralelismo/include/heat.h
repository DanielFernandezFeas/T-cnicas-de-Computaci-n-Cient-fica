//
// Created by DANIEL on 30/12/2025.
//

#ifndef T_CNICAS_DE_COMPUTACI_N_CIENT_FICA_HEAT_H
#define T_CNICAS_DE_COMPUTACI_N_CIENT_FICA_HEAT_H

#include "config.h"
#include "grid.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Perform a single heat diffusion step.
 * g: pointer to the current grid (updated in-place via swap)
 * g_copy: pointer to a temporary grid (may be NULL)
 * Returns 0 on success, negative on error.
 */
int heat_step(grid_t *g, grid_t *g_copy);

/*
 * Perform nsteps iterations of the heat solver on grid g.
 * Returns 0 on success, negative on error.
 */
int heat_iterate(grid_t *g, int nsteps);

#ifdef __cplusplus
}
#endif

#endif //T_CNICAS_DE_COMPUTACI_N_CIENT_FICA_HEAT_H