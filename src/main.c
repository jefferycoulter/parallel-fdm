#include "fdm.h"
#include "io.h"
#include "parallel.h"
#include "shape.h"

int main(int argc, char** argv)
{
    int n_proc, rank;
    int time_steps = 2000; // max time steps to iterate
    float radius = 20.0; // radius of region of interest

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialize the subdomain struct on each process
    Subdomain *sd = CreateSubdomain(n_proc, rank);

    // setup for finite difference computations
    DetermineStepSizes(sd);
    CreateShapeArray(sd, radius);
    SetBoundaryConditions(sd); 
    SetInitialConditions(sd);

    // computation
    ComputeFD(sd, Dirichlet, time_steps);

    SubdomainCleanUp(sd);
    MPI_Finalize();
    
    return 0;
}