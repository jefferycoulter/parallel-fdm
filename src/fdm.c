#include "fdm.h"
#include "io.h"
#include "shape.h"

#include <mpi.h>

#include <stdlib.h>
#include <math.h>
#include <memory.h>

//enum Mode{FD, Laplace};

void ComputeFD(Subdomain *sd, int bc, int time_steps)
{
    MPI_Gather(sd->u_now, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT, 
                   sd->u_global, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,
                   ROOT, MPI_COMM_WORLD);
    WriteData(*sd, FDM); // save initial conditions
    for (int t = 1; t < time_steps+1; t ++)
    {
        FTCS(sd, bc); // compute finite differences

        // send the local computations to ROOT process
        MPI_Gather(sd->u_next, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT, 
                   sd->u_global, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,
                   ROOT, MPI_COMM_WORLD);
        if (sd->rank == ROOT)
        {
            WriteData(*sd, FDM); // save each iteration
        }
    }
} // end void ComputeFD(Subdomain sd, int bc, int time_steps)

void FTCS(Subdomain *sd, int bc)
{
    int r = 0, c = 0, d = 0; // row, column, depth

    for (int i = 0; i < sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2]; i++)
    {
        InteriorOrBoundary(sd, bc, r, c, d)
        Increment(i, r, c, d)
    }

    // copy new data to old data array for next iteration
    // note:
    // sd.grid_l[2] = 1 for 2D case, so multiplying it here doesn't matter if the problem is 2D
    memcpy(sd->u_now, sd->u_next, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2] * sizeof((*sd->u_now)));
} // end void FTCS(Subdomain sd, int bc)

void InteriorFD(Subdomain *sd, int r, int c, int d)
{
    if (d == 0)
    {
        sd->u_next[r * sd->grid_l[1] + c] = FDX(sd, r, c, 0) + FDY(sd, r, c, 0);
    }
    else
    {
        sd->u_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = FDX(sd, r, c, d) + FDY(sd, r, c, d) + FDZ(sd, r, c, d);
    }
} // end void InteriorFD(Subdomain sd, int i, int j, int k)

void BoundaryFD(Subdomain *sd, int bc, int i, int j, int k)
{
    switch (bc)
    {
        case Dirichlet:
        {
            //SetBoundaryConditions(sd);
            // don't need to do anything assuming boundaries are constant
            break;
        } // end case Dirichlet
        case VonNeumann:
        {
            // need to implement
            break;
        } // end case VonNeumann
    } // end switch
} // end void BoundaryFD(Subdomain sd, int bc, int i, int j, int k)


void CreateShapeArray(Subdomain *subdomain, float radius)
{
    CoordShift(subdomain, radius);
    ApplyLaplaceFilter(subdomain);
} // end void CreateShapeArray(Subdomain subdomain, float radius)

void ApplyLaplaceFilter(Subdomain *sd)
{
    int r = 0, c = 0, d = 0; // row, column, depth

    for (int i = 0; i < sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2]; i++)
    {
        sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = Laplace(sd, r, c, d);
        AssignValue(sd, r, c, d);

        Increment(i, r, c, d)
    }
}

void SetBoundaryConditions(Subdomain *sd)
{
    int bc = 10.0;

    int r = 0, c = 0, d = 0;
    for (int i = 0; i < sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2]; i++)
    {
        switch (sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c])
        {
            case BOUNDARY:
                sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = bc;
                break;
            default:
                break;
        }

        Increment(i, r, c, d)
    }
} // end void SetBoundaryConditions(Subdomain sd)

void SetInitialConditions(Subdomain *sd)
{
    int ic = -10.0;

    int r = 0, c = 0, d = 0;
    for (int i = 0; i < sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2]; i++)
    {
        switch (sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c])
        {
            case INSIDE:
                sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = ic;
                break;
            default:
                break;
        }

        Increment(i, r, c, d)
    }
} // end void SetInitialConditions(Subdomain sd)

void DiscretizeSubdomain(Subdomain *sd)
{
    if(sd->grid_g[2] == 0)
    {
        sd->grid_g[0] = sd->dims_l[0] / sd->dx;
        sd->grid_g[1] = sd->dims_l[1] / sd->dy;
        sd->grid_g[2] = 0;
    }
    else
    {
        sd->grid_g[0] = sd->dims_l[0] / sd->dx;
        sd->grid_g[1] = sd->dims_l[1] / sd->dy;
        sd->grid_g[2] = sd->dims_l[2] / sd->dz;
    }
} // end void DiscretizeSubdomain(Subdomain sd)

void DetermineStepSizes(Subdomain *sd)
{
    // step sizes
    sd->dx = (float)sd->dims_g[0] / (float)sd->grid_g[0];
    sd->dy = (float)sd->dims_g[1] / (float)sd->grid_g[1];
    if(sd->n_dims == 3)
    {
        sd->dz = (float)sd->dims_g[2] / (float)sd->grid_g[2];
        sd->mu_z = sd->dt / pow(sd->dz, 2.0);
        if (sd->mu_z > 0.5) { fprintf(stderr, "Step sizes are too large. dz = %f, mu_z = %f > 0.5\n", sd->dx, sd->mu_z); exit(1); }
    }

    // stability factors
    sd->mu_x = sd->dt / pow(sd->dx, 2.0);
    sd->mu_y = sd->dt / pow(sd->dy, 2.0);


    // make sure step sizes are appropriate for stable results
    if (sd->mu_x > 0.5) { fprintf(stderr, "Step sizes are too large. dx = %f, mu_x = %f > 0.5\n", sd->dx, sd->mu_x); exit(1); }
    if (sd->mu_y > 0.5) { fprintf(stderr, "Step sizes are too large. dy = %f, mu_y = %f > 0.5\n", sd->dx, sd->mu_y); exit(1); }
} // end void DetermineStepSizes(Subdomain sd)


