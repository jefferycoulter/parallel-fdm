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
    Subdomain subdomain = { .n_proc = n_proc, .rank = rank,
                            .n_dims = 2, .dims_g = { 20, 20, 0 }, .grid_g = { 200, 200, 1 },
                            .dt = 0.05 // 50 milliseconds
                        };

    AllocateArraysFDM(&subdomain);
    PrepareSubdomains(&subdomain);
    CoordShift(&subdomain, radius);

    //MPI_Gather(subdomain.shape_arr_l, subdomain.grid_l[0] * subdomain.grid_l[1] * subdomain.grid_l[2], MPI_INT, 
    //               subdomain.shape_arr_g, subdomain.grid_l[0] * subdomain.grid_l[1] * subdomain.grid_l[2], MPI_INT,
    //               ROOT, MPI_COMM_WORLD);

    FILE *fp = fopen("data/shape.csv", "a");
    if (fp == NULL) { printf("file not found\n"); exit(1); }

    if (rank == ROOT)
    {
        //fprintf(stdout, "rank root\n");
        //fprintf(stdout, "grid_l[0] = %d\n", subdomain.grid_l[0]);
        //fprintf(stdout, "grid_l[1] = %d\n", subdomain.grid_l[1]);
        //fprintf(stdout, "grid_l[2] = %d\n", subdomain.grid_l[2]);
        for (int i = 0; i < subdomain.n_proc * subdomain.grid_l[0] * subdomain.grid_l[1] * subdomain.grid_l[2]; i++)
        {
            //fprintf(fp, "%d,", subdomain.shape_arr_g[i]);
        }
    }
    fprintf(fp, "\n");
    fclose(fp);
   
    free(subdomain.shape_arr_l);
    free(subdomain.shape_arr_g);
    free(subdomain.u_global);
    free(subdomain.u_now);
    free(subdomain.u_next);

    MPI_Finalize();

    return 0;
}