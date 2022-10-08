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
    //initialize the subdomain on each process
    Subdomain subdomain = { .n_proc = n_proc, .rank = rank,
                            .n_dims = 2, .dims_g = { 20, 20, 0},
                            .dx = 0.5, .dy = 0.5, .dz = 0.0, 
                            .dt = 0.05,
                            .mu_x = subdomain.dt / pow(subdomain.dx, 2.0),
                            .mu_y = subdomain.dt / pow(subdomain.dy, 2.0),
                            .mu_z = subdomain.dt / pow(subdomain.dz, 2.0)
                        };

    // get information related to subdomain position in the global context
    ComputeSubdomainDims(subdomain);
    GetSubdomainBounds(subdomain);

    // prepare subdomain for fdm
    DiscretizeSubdomain(subdomain);
    GenerateArraysFDM(subdomain);
    CreateShapeArray(subdomain, radius);

    SetBoundaryConditions(subdomain);
    SetInitialConditions(subdomain);

    ComputeFD(subdomain, Dirichlet, time_steps);

    free(subdomain.shape_arr);
    free(subdomain.u_now);
    free(subdomain.u_next);
    
    return 0;
}