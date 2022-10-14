#include "fdm.h"
#include "io.h"
#include "parallel.h"

#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int n_proc, rank;
    int time_steps = 5000; // max time steps to iterate
    float radius = 50.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialize the subdomain struct on each process
    Subdomain *sd = CreateSubdomain(n_proc, rank);

    // setup for finite difference computations
    DetermineStepSizes(sd);
    fprintf(stdout, "mu_x = %f\n", sd->mu_x);
    fprintf(stdout, "mu_y = %f\n", sd->mu_y);
    fprintf(stdout, "mu_z = %f\n", sd->mu_z);
    CreateShapeArray(sd, radius);
    
    SetBoundaryConditions(sd);
    SetInitialConditions(sd);

    //MPI_Gather(sd->u_now, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT, 
      //             sd->u_global, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,
        //           ROOT, MPI_COMM_WORLD);
    //WriteData((*sd), FDM);

    ComputeFD(sd, Dirichlet, time_steps);

    SubdomainCleanUp(sd);
   
    MPI_Finalize();
    
    return 0;
}