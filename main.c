#include "fdm.h"
#include "io.h"

#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    Data data = { 0 };
    // input corresponds to dx, dy, dz, dt
    StepSize step_sizes = { .dx = 0.5, .dy = 0.5, .dz = 0.0, .dt = 0.05 };
    Dimension dim = { .n = 2, .x = 20, .y = 20, .z = 0 };

    int time_steps = 1000;
    // number of finite difference cells or grid boxes along each dimension
    // pre-set to 3 as maximum number of dimensions and zeroed out so if 
    // no third dimension is used then it will be zero
    int cells[3] = { 0 }; 

    // generate subdomains from global dimension sizes before discretizing and computing

    DiscretizeSubdomain(dim, step_sizes, cells);
    GenerateSpatialArrays(&data, cells); // macro? also, should be done after creating subdomains

    SetBoundaryConditions(data, cells);
    SetInitialConditions(data.now, cells);

    ComputeFD(data, step_sizes, time_steps, cells);

    free(data.next);
    free(data.now);
    
    return 0;
}