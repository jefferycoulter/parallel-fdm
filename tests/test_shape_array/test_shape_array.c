#include "fdm.h"
#include "io.h"
#include "parallel.h"

#include <mpi.h>
#include <stdlib.h>

#define ROOT 0

int main(int argc, char **argv)
{
    int n_proc, rank;
    int time_steps = 1000; // max time steps to iterate
    float radius = 20.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //initialize the subdomain on each process
    Subdomain *subdomain = CreateSubdomain(n_proc, rank);

    CoordShift(subdomain, radius);

    MPI_Gather((*subdomain).shape_arr_l, (*subdomain).grid_l[0] * (*subdomain).grid_l[1] * (*subdomain).grid_l[2], MPI_INT, 
                   (*subdomain).shape_arr_g, (*subdomain).grid_l[0] * (*subdomain).grid_l[1] * (*subdomain).grid_l[2], MPI_INT,
                   ROOT, MPI_COMM_WORLD);

    if (rank == ROOT)
    {
        WriteData(*subdomain, Shape);
    }
   
    free((*subdomain).shape_arr_l);
    free((*subdomain).shape_arr_g);
    free((*subdomain).u_global);
    free((*subdomain).u_now);
    free((*subdomain).u_next);

    MPI_Finalize();

    return 0;
}