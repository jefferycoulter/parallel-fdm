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

    //CoordShift(sd, radius);
    fprintf(stdout, "rank %d: neighbor up: %d \t neighbor down: %d\n", sd->rank, sd->neighbors[0], sd->neighbors[1]);
    CreateShapeArray(sd, radius);

    CollectSubdomainData(sd);
    if ((*sd).rank == ROOT)
    {
        WriteData(*sd, Shape);
    }
   
    SubdomainCleanUp(sd);

    MPI_Finalize();

    return 0;
}