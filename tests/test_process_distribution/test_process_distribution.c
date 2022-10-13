#include "fdm.h"
#include "io.h"
#include "parallel.h"

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

    CoordShift(sd, radius);
    //CreateShapeArray(sd, radius);

    /*
    fprintf(stdout, "rank %d: sd->coords[0] = %d\n", sd->rank, sd->coords[0]);
    fprintf(stdout, "rank %d: sd->coords[1] = %d\n", sd->rank, sd->coords[1]);
    for (int n = 0; n < sd->n_dims; n++)
    {
        fprintf(stdout, "rank %d: sd->n_proc_dim[%d] = %d\n", sd->rank, n, sd->n_proc_dim2D[n]);
        fprintf(stdout, "rank %d: sd->grid_l[%d] = %d\n", sd->rank, n, sd->grid_l[n]);
        fprintf(stdout, "rank %d: sd->bounds_l[%d] = %d\n", sd->rank, 2*n, sd->bounds_l[2 * n]);
        fprintf(stdout, "rank %d: sd->bounds_l[%d] = %d\n", sd->rank, 2*n+1, sd->bounds_l[2 * n + 1]);
    }
    */
    CollectSubdomainData(sd);
    if ((*sd).rank == ROOT)
    {
        WriteData(*sd, Shape);
    }
   
    SubdomainCleanUp(sd);

    MPI_Finalize();

    return 0;
}