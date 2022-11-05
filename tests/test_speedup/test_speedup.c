#include "fdm.h"
#include "io.h"
#include "parallel.h"
#include "shape.h"

#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

int main(int argc, char** argv)
{
    struct timeval start_g, stop_g;
    struct timeval start_p, stop_p;
    double time_p, time_g;

    int n_proc, rank;
    int time_steps = 4000; // max time steps to iterate
    float radius = 45.0; // radius of region of interest

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == ROOT) { gettimeofday(&start_g, NULL); }

    // initialize the subdomain struct on each process
    Subdomain *sd = CreateSubdomain(n_proc, rank);

    // setup for finite difference computations
    DetermineStepSizes(sd);
    if (rank == ROOT) { gettimeofday(&start_p, NULL); }
    CreateShapeArray(sd, radius); // parallel
    SetBoundaryConditions(sd); // parallel
    SetInitialConditions(sd); // parallel

    // computation
    ComputeFD(sd, Dirichlet, time_steps); // parallel

    if (rank == ROOT) { gettimeofday(&stop_p, NULL); }
    if (rank == ROOT)
    { 
        //time_p = stop_p.tv_sec - start_p.tv_sec;
        //FILE *fp = fopen("data/times.csv", "a");
        //if (fp == NULL) { printf("file not found\n"); exit(1); }
        //fprintf(fp, "%d,", sd->n_proc);
        //fprintf(fp, "%f,", time_p);
        //fprintf(fp, "\n");
        //fclose(fp);
    }

    SubdomainCleanUp(sd);

    if (rank == ROOT) { gettimeofday(&stop_g, NULL); }

    if (rank == ROOT)
    {
        time_g = stop_g.tv_sec - start_g.tv_sec;
        time_p = stop_p.tv_sec - start_p.tv_sec;
        fprintf(stdout, "global time: %f\n", time_g);
        fprintf(stdout, "parallelizable time: %f\n", time_p);
    }

    MPI_Finalize();
    
    return 0;
}