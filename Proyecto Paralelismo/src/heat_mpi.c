#include <mpi.h>
#include <stdlib.h>
#include "heat.h"
#include "grid.h"

int heat_step_mpi(grid_t *local,
                  grid_t *tmp,
                  int local_ny,
                  MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int nx = local->nx;

    //intercambia filas de halo con vecinos
    
    if (rank > 0)
    {
        MPI_Sendrecv(
            &GRID_AT(local, 1, 0), nx, MPI_REAL,
            rank - 1, 0,
            &GRID_AT(local, 0, 0), nx, MPI_REAL,
            rank - 1, 0,
            comm, MPI_STATUS_IGNORE);
    }

    if (rank < size - 1)
    {
        MPI_Sendrecv(
            &GRID_AT(local, local_ny, 0), nx, MPI_REAL,
            rank + 1, 0,
            &GRID_AT(local, local_ny + 1, 0), nx, MPI_REAL,
            rank + 1, 0,
            comm, MPI_STATUS_IGNORE);
    }

    // actualiza la malla local (no tocar las filas de contorno globales)
    for (int y = 1; y <= local_ny; ++y)
    {
        // Skip global top boundary (owned by rank 0) and
        // global bottom boundary (owned by rank size-1).
        if ((rank == 0 && y == 1) || (rank == size - 1 && y == local_ny))
            continue;

        for (int x = 1; x < nx - 1; ++x)
        {
            GRID_AT(tmp, y, x) =
                (GRID_AT(local, y - 1, x) +
                 GRID_AT(local, y + 1, x) +
                 GRID_AT(local, y, x - 1) +
                 GRID_AT(local, y, x + 1)) *
                (real_t)0.25;
        }
    }

    grid_swap(local, tmp);
    return 0;
}

/* ------------------------------------------------------------
 * Iteraciones MPI (global)
 * ------------------------------------------------------------ */
int heat_iterate_mpi(grid_t *global,
                     int nsteps,
                     MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int nx = global->nx;
    int ny = global->ny;

    int base = ny / size;
    int rem = ny % size;

    int local_ny = base + (rank < rem ? 1 : 0);

    // reserva de mallas locales y halos
    grid_t local = {0}, tmp = {0};
    if (grid_alloc(&local, nx, local_ny + 2) != 0)
        return -2;
    if (grid_alloc(&tmp, nx, local_ny + 2) != 0)
        return -2;

    // reparto de filas entre procesos
    int *sendcounts = NULL, *displs = NULL;

    if (rank == 0)
    {
        sendcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));

        int off = 0;
        for (int r = 0; r < size; ++r)
        {
            int rows = base + (r < rem);
            sendcounts[r] = rows * nx;
            displs[r] = off;
            off += rows * nx;
        }
    }

    // reparto inicial con Scatterv
    MPI_Scatterv(global->data, sendcounts, displs, MPI_REAL,
                 &GRID_AT(&local, 1, 0),
                 local_ny * nx, MPI_REAL,
                 0, comm);

    for (int it = 0; it < nsteps; ++it)
    {
        heat_step_mpi(&local, &tmp, local_ny, comm);
    }

    // recolección final con Gatherv
    MPI_Gatherv(&GRID_AT(&local, 1, 0),
                local_ny * nx, MPI_REAL,
                global->data, sendcounts, displs, MPI_REAL,
                0, comm);

    if (rank == 0)
    {
        free(sendcounts);
        free(displs);
    }

    grid_free(&local);
    grid_free(&tmp);
    return 0;
}
