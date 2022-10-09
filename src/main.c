#include "fdm.h"
#include "io.h"
#include "parallel.h"

#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n_proc, rank;
    int time_steps = 1000; // max time steps to iterate
    float radius;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialize the subdomain struct on each process
    Subdomain subdomain = { .n_proc = n_proc, .rank = rank,
                            .n_dims = 2, .dims_g = { 20, 20, 0 }, .grid_g = { 100, 100, 1 },
                            .dt = 0.05 // 50 milliseconds
                        };

    // perform further initializations and allocations
    PrepareSubdomains(&subdomain);

    // setup for finite difference computations
    DetermineStepSizes(subdomain);
    AllocateArraysFDM(&subdomain);
    CreateShapeArray(subdomain, radius);

    SetBoundaryConditions(subdomain);
    SetInitialConditions(subdomain);

    ComputeFD(subdomain, Dirichlet, time_steps);
   
    MPI_Finalize();
    
    return 0;
}