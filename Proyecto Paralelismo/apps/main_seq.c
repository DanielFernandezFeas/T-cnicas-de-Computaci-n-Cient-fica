#include <stdio.h>
#include "grid.h"
#include "config.h"


int main(void)
{
    grid_t g;
    if (grid_alloc(&g, 10, 10) != 0)
    {
        printf("Error reservando malla\n");
        return 1;
    }

    grid_init_bc(&g, 100.0, 0.0);
    printf("Checksum: %.2f\n", (double)grid_checksum(&g));
    grid_print(&g, 10, 10);

    grid_free(&g);
    return 0;
}
