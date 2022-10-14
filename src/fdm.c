#include "fdm.h"
#include "io.h"

#include <mpi.h>

#include <stdlib.h>
#include <math.h>
#include <memory.h>

enum Mode{FD, Laplace};

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
    for (int i = 0; i < sd->grid_l[0]; i ++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            if (sd->n_dims == 2) // 2D problem
            {   
                switch (sd->shape_next[i * sd->grid_l[1] + j])
                {
                case OUTSIDE:
                    break;
                case INSIDE:
                    InteriorFD(sd, i, j, 0);
                    break;
                case BOUNDARY:
                    BoundaryFD(sd, bc, i, j, 0);
                    break;
                }
            }
            else if (sd->n_dims == 3)// 3D problem
            {
                for (int k = 0; k < sd->grid_l[2]; k++)
                {
                    switch (sd->shape_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k])
                    {
                    case OUTSIDE:
                        break;
                    case INSIDE:
                        InteriorFD(sd, i, j, k);
                        break;
                    case BOUNDARY:
                        BoundaryFD(sd, bc, i, j, k);
                        break;
                    } // end switch
                } // end k loop
            }
            else
            {
                fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
                exit(1);
            } // end if-else
        } // end j loop
    } // end i loop

    // copy new data to old data array for next iteration
    // note:
    // sd.grid_l[2] = 1 for 2D case, so multiplying it here doesn't matter if the problem is 2D
    memcpy(sd->u_now, sd->u_next, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2] * sizeof((*sd->u_now)));
} // end void FTCS(Subdomain sd, int bc)

void InteriorFD(Subdomain *sd, int i, int j, int k)
{
    if (k == 0) // 2D problem
    {
        sd->u_next[i * sd->grid_l[1] + j] = (1 - 2 * sd->mu_x - 2 * sd->mu_y) * sd->u_now[i * sd->grid_l[1] + j]          \
            + sd->mu_x * (sd->u_now[i * sd->grid_l[1] + j + sd->grid_l[0]] + sd->u_now[i * sd->grid_l[1] + j - sd->grid_l[0]])  \
            + sd->mu_y * (sd->u_now[i * sd->grid_l[1] + j + 1] + sd->u_now[i * sd->grid_l[1] + j - 1]);
    }
    else
    {

            sd->u_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = (1 - 2 * sd->mu_x - 2 * sd->mu_y - 2 * sd->mu_z) * sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k]                   \
            + sd->mu_x * (sd->u_now[(i * sd->grid_l[1] + j + sd->grid_l[0]) * sd->grid_l[2] + k] + sd->u_now[(i * sd->grid_l[1] + j - sd->grid_l[0]) * sd->grid_l[2] + k])    \
            + sd->mu_y * (sd->u_now[(i * sd->grid_l[1] + j + 1) * sd->grid_l[2] + k] + sd->u_now[(i * sd->grid_l[1] + j - 1) * sd->grid_l[2] + k])                  \
            + sd->mu_z * (sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k + sd->grid_l[0] * sd->grid_l[1]] + sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k - sd->grid_l[0] * sd->grid_l[1]]);
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
            if (k == 0) // 2D problem
            {
                sd->u_next[i * sd->grid_l[1] + j] = (1 - 2 * sd->mu_x - 2 * sd->mu_y) * sd->u_now[i * sd->grid_l[1] + j]          \
                    + sd->mu_x * (sd->u_now[i * sd->grid_l[1] + j + sd->grid_l[0]] + sd->u_now[i * sd->grid_l[1] + j - sd->grid_l[0]])  \
                    + sd->mu_y * (sd->u_now[i * sd->grid_l[1] + j + 1] + sd->u_now[i * sd->grid_l[1] + j - 1]);
            }
            else // 3d problem
            {
                sd->u_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = (1 - 2 * sd->mu_x - 2 * sd->mu_y - 2 * sd->mu_z) * sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k]                   \
                + sd->mu_x * (sd->u_now[(i * sd->grid_l[1] + j + sd->grid_l[0]) * sd->grid_l[2] + k] + sd->u_now[(i * sd->grid_l[1] + j - sd->grid_l[0]) * sd->grid_l[2] + k])    \
                + sd->mu_y * (sd->u_now[(i * sd->grid_l[1] + j + 1) * sd->grid_l[2] + k] + sd->u_now[(i * sd->grid_l[1] + j - 1) * sd->grid_l[2] + k])                  \
                + sd->mu_z * (sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k + sd->grid_l[0] * sd->grid_l[1]] + sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k - sd->grid_l[0] * sd->grid_l[1]]);
            }
            if (sd->grid_g[2] == 0)
            {
                // corners
                // bottom left corner
                // bottom right corner
                // top left corner
                // top right corner

                // edges
                // left edge
                // right edge
                // bottom edge
                // top edge
            }
            else
            {
                    // corners
                    // front bottom left corner
                    // front bottom right corner
                    // front top left corner
                    // front top right corner
                    // back bottom left corner
                    // back bottom right corner
                    // back top left corner
                    // back top right corner

                    // edges
                    // front left edge
                    // front right edge
                    // front bottom edge
                    // front top edge
                    // back left edge
                    // back right edge
                    // back bottom edge
                    // back top edge
                    // middle left edge
                    // middle right edge
                    // middle bottom edge
                    // middle top edge
            } // end if-else
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
    for (int i = 1; i < sd->grid_l[0] - 1; i++)
    {
        for (int j = 1; j < sd->grid_l[1] - 1; j++)
        {
            if (sd->n_dims == 2)
            {
                // apply laplace filter to each point
                sd->shape_next[i * sd->grid_l[1] + j] = -4 * sd->shape_now[i * sd->grid_l[1] + j]          \
                    + sd->shape_now[i * sd->grid_l[1] + j + sd->grid_l[0]] + sd->shape_now[i * sd->grid_l[1] + j - sd->grid_l[0]]  \
                    + sd->shape_now[i * sd->grid_l[1] + j + 1] + sd->shape_now[i * sd->grid_l[1] + j - 1];

                // assign point correct value for finite difference computations
                switch (sd->shape_next[i * sd->grid_l[1] + j])
                {
                    case 0: // if the new value is zero, then the previous value was either INSIDE or OUTSIDE
                        switch (sd->shape_now[i * sd->grid_l[1] + j])
                        {
                            case OUTSIDE: // if previous value was OUTSIDE, then new value is still outside
                                sd->shape_next[i * sd->grid_l[1] + j] = OUTSIDE;
                                break;
                            
                            case INSIDE: // if previous value was INSIDE, then switch it back to inside
                                sd->shape_next[i * sd->grid_l[1] + j] = INSIDE;
                                break;
                        } // end switch
                        break;
                    default: // if the new value is not zero, then it is a boundary
                        sd->shape_next[i * sd->grid_l[1] + j] = BOUNDARY;
                        break;
                } // end switch
            }
            else if (sd->n_dims == 3)
            {
                for (int k = 1; k < sd->grid_l[2] - 1; k ++)
                {
                    // apply laplace filter to each point
                    sd->shape_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = -6 * sd->shape_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k]                   \
                    + sd->shape_now[(i * sd->grid_l[1] + j + sd->grid_l[0]) * sd->grid_l[2] + k] + sd->shape_now[(i * sd->grid_l[1] + j - sd->grid_l[0]) * sd->grid_l[2] + k]   \
                    + sd->shape_now[(i * sd->grid_l[1] + j + 1) * sd->grid_l[2] + k] + sd->shape_now[(i * sd->grid_l[1] + j - 1) * sd->grid_l[2] + k]                  \
                    + sd->shape_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k + sd->grid_l[0] * sd->grid_l[1]] + sd->shape_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k - sd->grid_l[0] * sd->grid_l[1]];

                    // assign point correct value for finite difference computations
                    switch (sd->shape_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k])
                    {
                        case 0: // if the new value is zero, then the previous value was either INSIDE or OUTSIDE
                            switch (sd->shape_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k])
                            {
                                case OUTSIDE: // if previous value was OUTSIDE, then new value is still outside
                                    sd->shape_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = OUTSIDE;
                                    break;
                                
                                case INSIDE: // if previous value was INSIDE, then switch it back to inside
                                    sd->shape_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = INSIDE;
                                    break;
                            }
                            break;
                        default: // if the new value is not zero, then it is a boundary
                            sd->shape_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = BOUNDARY;
                            break;
                    }
                }
            }
            else
            {
                fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
                exit(1);
            }
        }
    }
}

void SetBoundaryConditions(Subdomain *sd)
{
    int bc = 10.0;

    for (int i = 0; i < sd->grid_l[0]; i++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            if (sd->n_dims == 2)
            {
                switch (sd->shape_next[i * sd->grid_l[1] + j])
                {
                    case BOUNDARY:
                        sd->u_now[i * sd->grid_l[1] + j] = bc;
                        break;
                    default:
                        break;
                }
            }
            else if (sd->n_dims == 3)
            {
                for (int k = 0; k < sd->grid_l[2]; k++)
                {

                }
            }
            else
            {
                fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
                exit(1);
            }
        }
    }

} // end void SetBoundaryConditions(Subdomain sd)

void SetInitialConditions(Subdomain *sd)
{
    int ic = -10.0;

    for (int i = 0; i < sd->grid_l[0]; i++)
    {
        for (int j = 0; j < sd->grid_l[1]; j++)
        {
            if (sd->n_dims == 2)
            {
                switch (sd->shape_next[i * sd->grid_l[1] + j])
                {
                    case INSIDE:
                        sd->u_now[i * sd->grid_l[1] + j] = ic;
                        break;
                    default:
                        break;
                }
            }
            else if (sd->n_dims == 3)
            {
                for (int k = 0; k < sd->grid_l[2]; k++)
                {

                }
            }
            else
            {
                fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
                exit(1);
            }
        }
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


