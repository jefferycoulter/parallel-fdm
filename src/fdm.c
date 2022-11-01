#include "fdm.h"
#include "io.h"
#include "shape.h"

#include <stdlib.h>
#include <math.h>
#include <memory.h>

void ComputeFD(Subdomain *sd, int bc, int time_steps)
{
    CollectSubdomainData(sd, FDM, 0);
    if (sd->rank == ROOT) { WriteData(*sd, FDM); }
    ShareGhosts(sd, FDM);
    for (int time = 1; time < time_steps+1; time ++)
    {
        FTCS(sd, bc); // forward-time central-space scheme
        ShareGhosts(sd, FDM);
        if (time % 50 == 0)
        {
            CollectSubdomainData(sd, FDM, time);
            if (sd->rank == ROOT) { WriteData(*sd, FDM); }
        }
    }
} // end void ComputeFD(Subdomain *sd, int bc, int time_steps)

void FTCS(Subdomain *sd, int bc)
{
    int id;
    int offset = sd->ghost_size;

    for (int i = 0; i < sd->grid_l[0]; i++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            for (int k = 0; k < sd->grid_l[2]; k++) // for 2d case grid_l[2] = 1, k will always be zero
            {
                id = offset + ID(sd, i, j, k);
                switch (sd->shape_next[id])
                    {
                    case Outside:
                        break;
                    case Inside:
                        InteriorFD(sd, id);
                        break;
                    case Boundary:
                        BoundaryFD(sd, bc, id);
                        break;
                    } // end switch
            } // end k loop
        } // end j loop
    } // end i loop

    // copy new data to old data array for next iteration
    // note:
    // sd.grid_l[2] = 1 for 2D case, so multiplying it here doesn't matter if the problem is 2D
    memcpy(sd->u_now, sd->u_next, sizeof(float) * ((*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] + (2 * (*sd).ghost_size)));
} // end void FTCS(Subdomain *sd, int bc)

void InteriorFD(Subdomain *sd, int id)
{
    sd->u_next[id] = FD(sd, id);
} // end void InteriorFD(Subdomain *sd, int id)

void BoundaryFD(Subdomain *sd, int bc, int id)
{
    switch (bc)
    {
        case Dirichlet:
        {
            sd->u_next[id] = sd->u_now[id];
            break;
        } // end case Dirichlet
        case VonNeumann:
        {
            // need to implement
            break;
        } // end case VonNeumann
    } // end switch
} // end void BoundaryFD(Subdomain *sd, int bc, int id)

void CreateShapeArray(Subdomain *sd, float radius)
{
    CoordShift(sd, radius);

    ShareGhosts(sd, Shape);
    MPI_Barrier(sd->COMM_FDM);

    ApplyLaplaceFilter(sd);
    ShareGhosts(sd, Shape);
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
                sd->shape_next[id] = Laplace(sd, id);
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
                id = offset + ID(sd, i, j, k);
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
                id = offset + ID(sd, i, j, k);
                switch (sd->shape_next[id])
                {
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

void DetermineStepSizes(Subdomain *sd)
{
    float D = 3; // diffusion constant
    // step sizes
    sd->dx = (float)sd->dims_g[0] / (float)sd->grid_g[0];
    sd->dy = (float)sd->dims_g[1] / (float)sd->grid_g[1];
    sd->dz = (float)sd->dims_g[2] / (float)sd->grid_g[2];

    // stability factors
    sd->mu_x = (sd->dt / pow(sd->dx, 2.0)) * D;
    sd->mu_y = (sd->dt / pow(sd->dy, 2.0)) * D;
    sd->mu_z = (sd->dt / pow(sd->dz, 2.0)) * D;

    // make sure step sizes are appropriate for stable results
    if (sd->mu_x > 0.125) { fprintf(stderr, "Step sizes are too large. dx = %.3f, mu_x = %.3f > 0.125\n", sd->dx, sd->mu_x); exit(1); }
    if (sd->mu_y > 0.125) { fprintf(stderr, "Step sizes are too large. dy = %.3f, mu_y = %.3f > 0.125\n", sd->dx, sd->mu_y); exit(1); }
    if (sd->mu_z > 0.125) { fprintf(stderr, "Step sizes are too large. dz = %.3f, mu_z = %.3f > 0.125\n", sd->dx, sd->mu_z); exit(1); }
} // end void DetermineStepSizes(Subdomain sd)


