#include "config.h"
#include "grid.h"

#ifdef USE_MPI
  #include <mpi.h>
  int heat_step_mpi(grid_t *g, grid_t *g_copy, int local_ny, MPI_Comm comm);
int heat_iterate_mpi(grid_t *g, int nsteps, MPI_Comm comm);
#endif

int heat_step_seq(grid_t *g, grid_t *g_copy);
int heat_iterate_seq(grid_t *g, int nsteps);
int heat_step_omp(grid_t *g, grid_t *g_copy, int nthreads);
int heat_iterate_omp(grid_t *g, int nsteps, int nthreads);
