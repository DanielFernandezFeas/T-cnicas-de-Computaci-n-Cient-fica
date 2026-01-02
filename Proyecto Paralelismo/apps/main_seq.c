#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "config.h"
#include "heat.h"


int main(int argc, char **argv)
{

    int nx, ny, iters;
    grid_t g;

    if (argc > 2)
    {
        nx = atoi(argv[1]);
        ny = atoi(argv[2]);
    }
    else
    {
        nx = 100;
        ny = 100;
    }

    if(argc > 3)
    {
        iters = atoi(argv[3]);
    }
    else
    {
        iters = 1000;
    }


    if (grid_alloc(&g, nx, ny) != 0)
    {
        printf("Error reservando malla\n");
        return 1;
    }


    grid_init_bc(&g, 100.0, 0.0);
    printf("Malla inicial %dx%d:\n", nx, ny);
    printf("Checksum: %.2f\n", (double)grid_checksum(&g));
    grid_print(&g, 10, 10);
    grid_apply_bc(&g, 150.0, 25.0);
    heat_iterate(&g, iters);
    printf("Checksum despues de iterar: %.2f\n", (double)grid_checksum(&g));
    grid_print(&g, 10, 10);
    grid_free(&g);
    return 0;
}
