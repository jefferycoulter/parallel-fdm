#include "fdm.h"
#include "io.h"
#include "parallel.h"
#include "shape.h"

#include <mpi.h>
#include <stdlib.h>

int main(int argc, char **argv) 
{
    int n_proc, rank;
    MPI_Status status;
    int time_steps = 1000; // max time steps to iterate
    float radius = 40.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //initialize the subdomain on each process
    Subdomain *sd = CreateSubdomain(n_proc, rank);

    // stability factors for finite difference
    Mu mu = {   .x = (float)sd->dims_g[0] / (float)sd->grid_g[0],
                .y = (float)sd->dims_g[1] / (float)sd->grid_g[1],
                .z = (float)sd->dims_g[2] / (float)sd->grid_g[2],
    };

    Mu mu_s = { .x = 1, .y = 1, .z = 1 };
    
    int (*shape_now)[(*sd).grid_l[1]][(*sd).grid_l[2]] = malloc(sizeof(int[(*sd).grid_l[0]][(*sd).grid_l[1]][(*sd).grid_l[2]]));
    int (*shape_next)[(*sd).grid_l[1]][(*sd).grid_l[2]] = malloc(sizeof(int[(*sd).grid_l[0]][(*sd).grid_l[1]][(*sd).grid_l[2]]));
    float (*u_now)[(*sd).grid_l[1]][(*sd).grid_l[2]] = malloc(sizeof(float[(*sd).grid_l[0]][(*sd).grid_l[1]][(*sd).grid_l[2]]));
    float (*u_next)[(*sd).grid_l[1]][(*sd).grid_l[2]] = malloc(sizeof(float[(*sd).grid_l[0]][(*sd).grid_l[1]][(*sd).grid_l[2]]));
    
    int (*shape_g)[(*sd).grid_g[1]][(*sd).grid_g[2]]  = malloc(sizeof(int[(*sd).grid_g[0]][(*sd).grid_g[1]][(*sd).grid_g[2]]));
    float *u_g; // global solution array

    CreateShapeArray(sd, shape_next, shape_now, mu_s, radius);

    CollectSubdomainData(sd, shape_next, shape_g);
    if ((*sd).rank == ROOT)
    {
        WriteData(sd, shape_g, Shape);
    }
   
    //SubdomainCleanUp(sd);

    MPI_Finalize();

    return 0;
}