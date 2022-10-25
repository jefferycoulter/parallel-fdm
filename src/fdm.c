#include "fdm.h"
#include "io.h"
#include "shape.h"

#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include <mpi.h>

void ComputeFD(Subdomain *sd, int bc, int time_steps)
{
    MPI_Gather(sd->u_now, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT, 
                   sd->u_global, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT,
                   ROOT, MPI_COMM_WORLD);
    int t = 0;
    WriteData(*sd, FDM); // save initial conditions
    for (int t = 1; t < time_steps+1; t ++)
    {
        FTCS(sd, bc); // forward-time central-space scheme

        // send the local computations to ROOT process
        MPI_Gather(sd->u_next, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT, 
                   sd->u_global, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT,
                   ROOT, MPI_COMM_WORLD);
        if (sd->rank == ROOT)
        {
            WriteData(*sd, FDM); // save each iteration
        }
    }
} // end void ComputeFD(Subdomain *sd, int bc, int time_steps)

void FTCS(Subdomain *sd, int bc)
{
    int id;
    int offset = sd->ghost_size;

    for (int i = 0; i < sd->grid_l[0]; i ++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            for (int k = 0; k < sd->grid_l[2]; k++) // for 2d case grid_l[2] = 1, k will always be zero
            {
                id = ID(sd, i, j, k) + offset;
                switch (sd->shape_next[id])
                    {
                    case Outside:
                        break;
                    case Inside:
                        InteriorFD(sd, i, j, k);
                        break;
                    case Boundary:
                        BoundaryFD(sd, bc, i, j, k);
                        break;
                    } // end switch
            } // end k loop
        } // end j loop
    } // end i loop

    // copy new data to old data array for next iteration
    // note:
    // sd.grid_l[2] = 1 for 2D case, so multiplying it here doesn't matter if the problem is 2D
    memcpy(sd->u_now, sd->u_next, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof((*sd->u_now)));
} // end void FTCS(Subdomain *sd, int bc)

void InteriorFD(Subdomain *sd, int i, int j, int k)
{
    int offset = sd->ghost_size;
    int id = ID(sd, i, j, k) + offset;
    sd->u_next[id] = FD(sd, i, j, k, offset);
} // end void InteriorFD(Subdomain *sd, int i, int j, int k)

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
} // end void BoundaryFD(Subdomain *sd, int bc, int i, int j, int k)


void CreateShapeArray(Subdomain *sd, float radius)
{
    CoordShift(sd, radius);
    ShareGhosts(sd, Shape);

    MPI_Barrier(sd->COMM_FDM);
    ApplyLaplaceFilter(sd);
} // end void CreateShapeArray(Subdomain *sd, float radius)

void ApplyLaplaceFilter(Subdomain *sd)
{
    int id; // index in linear memory
    int offset = sd->ghost_size; // offset if running parallel

    for (int i = 0; i < sd->grid_l[0]; i++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            for (int k = 0; k < sd->grid_l[2]; k++)
            {
                id = offset + ID(sd, i, j, k);
                sd->shape_next[id] = Laplace(sd, i, j, k, offset);
                AssignValue(sd, id)
            } // end k loop
        } // end j loop
    } // end i loop
} // end void ApplyLaplaceFilter(Subdomain *sd)

void SetBoundaryConditions(Subdomain *sd)
{
    int id;
    int offset = sd->ghost_size;

    int bc = 10.0;

    for (int i = 0; i < sd->grid_l[0]; i++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            for (int k = 0; k < sd->grid_l[2]; k ++)
            {
                id = ID(sd, i, j, k) + offset;
                switch (sd->shape_next[id])
                {
                    case Boundary:
                        sd->u_now[id] = bc;
                        break;
                    default:
                        break;
                } // end switch
            } // end k loop
        } // end j loop
    } // end i loop
} // end void SetBoundaryConditions(Subdomain sd)

void SetInitialConditions(Subdomain *sd)
{
    int id;
    int offset = sd->ghost_size;

    int ic = -10.0;

    for (int i = 0; i < sd->grid_l[0]; i++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            for (int k = 0; k < sd->grid_l[2]; k++)
            {
                switch (sd->shape_next[id])
                {
                    id = ID(sd, i, j, k) + offset;
                    case Inside:
                        sd->u_now[id] = ic;
                        break;
                    default:
                        break;
                } // end switch
            } // end k loop
        } // end j loop
    } // end i loop
} // end void SetInitialConditions(Subdomain sd)

void DiscretizeSubdomain(Subdomain *sd)
{
    if(sd->grid_g[2] == 0)
    {
        sd->grid_g[0] = sd->dims_g[0] / sd->dx;
        sd->grid_g[1] = sd->dims_g[1] / sd->dy;
        sd->grid_g[2] = 0;
    }
    else
    {
        sd->grid_g[0] = sd->dims_g[0] / sd->dx;
        sd->grid_g[1] = sd->dims_g[1] / sd->dy;
        sd->grid_g[2] = sd->dims_g[2] / sd->dz;
    }
} // end void DiscretizeSubdomain(Subdomain sd)

void DetermineStepSizes(Subdomain *sd)
{
    float D = 5; // diffusion constant
    // step sizes
    sd->dx = (float)sd->dims_g[0] / (float)sd->grid_g[0];
    sd->dy = (float)sd->dims_g[1] / (float)sd->grid_g[1];
    sd->dz = (float)sd->dims_g[2] / (float)sd->grid_g[2];

    // stability factors
    sd->mu_x = (sd->dt / pow(sd->dx, 2.0)) * D;
    sd->mu_y = (sd->dt / pow(sd->dy, 2.0)) * D;
    sd->mu_z = (sd->dt / pow(sd->dz, 2.0)) * D;

    // make sure step sizes are appropriate for stable results
    if (sd->mu_x > 0.5) { fprintf(stderr, "Step sizes are too large. dx = %f, mu_x = %f > 0.5\n", sd->dx, sd->mu_x); exit(1); }
    if (sd->mu_y > 0.5) { fprintf(stderr, "Step sizes are too large. dy = %f, mu_y = %f > 0.5\n", sd->dx, sd->mu_y); exit(1); }
    if (sd->mu_z > 0.5) { fprintf(stderr, "Step sizes are too large. dz = %f, mu_z = %f > 0.5\n", sd->dx, sd->mu_z); exit(1); }
} // end void DetermineStepSizes(Subdomain sd)


