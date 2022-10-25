#include "fdm.h"
#include "io.h"
#include "parallel.h"
#include "shape.h"

#include <mpi.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    int n_proc, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //initialize the subdomain on each process
    Subdomain *sd = CreateSubdomain(n_proc, rank);
    
    fprintf(stdout, "rank %d: coods = (%d, %d, %d)\n",  sd->rank,       \
                                                        sd->coords[0],  \
                                                        sd->coords[1],  \
                                                        sd->coords[2]);

    fprintf(stdout, "rank %d: n_proc_dim = (%d, %d, %d)\n", sd->rank,           \
                                                            sd->n_proc_dim[0],  \
                                                            sd->n_proc_dim[1],  \
                                                            sd->n_proc_dim[2]);
    fprintf(stdout, "rank %d: bounds = [%d:%d, %d:%d, %d:%d]\n",  sd->rank,                         \
                                                                sd->bounds_l[0], sd->bounds_l[1],   \
                                                                sd->bounds_l[2], sd->bounds_l[3],   \
                                                                sd->bounds_l[4], sd->bounds_l[5]);
    fprintf(stdout, "rank %d: neighbor up: %d \t neighbor down: %d\n",  sd->rank,               \
                                                                        sd->neighbors[0],       \
                                                                        sd->neighbors[1]);
   
    SubdomainCleanUp(sd);

    MPI_Finalize();

    return 0;
}