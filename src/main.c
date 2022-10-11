#include "fdm.h"
#include "io.h"
#include "parallel.h"

#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n_proc, rank;
    MPI_Status status;
    int time_steps = 1000; // max time steps to iterate
    float radius;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialize the subdomain struct on each process
    Subdomain *subdomain = CreateSubdomain(n_proc, rank);

    // perform further initializations and allocations

    // setup for finite difference computations
    DetermineStepSizes(subdomain);
    CreateShapeArray(subdomain, radius);

    SetBoundaryConditions(subdomain);
    SetInitialConditions(subdomain);

    ComputeFD(subdomain, Dirichlet, time_steps);

   
    MPI_Finalize();
    
    return 0;
}