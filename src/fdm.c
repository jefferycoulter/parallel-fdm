#include "fdm.h"
#include "io.h"

#include <mpi.h>

#include <stdlib.h>
#include <math.h>
#include <memory.h>

void ComputeFD(Subdomain *sd, int bc, int time_steps)
{
    WriteData(*sd, FDM); // save initial conditions
    for (int t = 1; t < time_steps+1; t ++)
    {
        FTCS(sd, bc); // compute difference

        // send the local computations to ROOT process
        MPI_Gather(sd->u_next, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT, 
                   sd->u_global, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,
                   ROOT, MPI_COMM_WORLD);
        if (sd->rank == ROOT)
        {
            WriteData(*sd, FDM); // save each iteration
        }
    }

    // free memory once finished
    free((*sd).shape_arr_l);
    free((*sd).shape_arr_g);
    free((*sd).u_global);
    free((*sd).u_now);
    free((*sd).u_next);
} // end void ComputeFD(Subdomain sd, int bc, int time_steps)

void FTCS(Subdomain *sd, int bc)
{
    for (int i = 0; i < sd->grid_l[0]; i ++)
    {
        for (int j = 0; j < sd->grid_l[0]; j++)
        {
            if (sd->grid_l[2] == 0) // 2D problem
            {
                switch (sd->shape_arr_l[i * sd->grid_l[1] + j])
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
            else // 3D problem
            {
                for (int k = 0; k < sd->grid_l[2]; k++)
                {
                    switch (sd->shape_arr_l[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k])
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
            // don't need to do anything assuming boundaries are constant
            break;
        } // end case Dirichlet
        case VonNeumann:
        {
            if (k == 0) // 2D problem
            {
                sd->u_next[i * sd->grid_g[1] + j] = (1 - 2 * sd->mu_x - 2 * sd->mu_y) * sd->u_now[i * sd->grid_g[1] + j]          \
                    + sd->mu_x * (sd->u_now[i * sd->grid_g[1] + j + sd->grid_g[0]] + sd->u_now[i * sd->grid_g[1] + j - sd->grid_g[0]])  \
                    + sd->mu_y * (sd->u_now[i * sd->grid_g[1] + j + 1] + sd->u_now[i * sd->grid_g[1] + j - 1]);
            }
            else // 3d problem
            {
                sd->u_next[(i * sd->grid_g[1] + j) * sd->grid_g[2] + k] = (1 - 2 * sd->mu_x - 2 * sd->mu_y - 2 * sd->mu_z) * sd->u_now[(i * sd->grid_g[1] + j) * sd->grid_g[2] + k]                   \
                + sd->mu_x * (sd->u_now[(i * sd->grid_g[1] + j + sd->grid_g[0]) * sd->grid_g[2] + k] + sd->u_now[(i * sd->grid_g[1] + j - sd->grid_g[0]) * sd->grid_g[2] + k])    \
                + sd->mu_y * (sd->u_now[(i * sd->grid_g[1] + j + 1) * sd->grid_g[2] + k] + sd->u_now[(i * sd->grid_g[1] + j - 1) * sd->grid_g[2] + k])                  \
                + sd->mu_z * (sd->u_now[(i * sd->grid_g[1] + j) * sd->grid_g[2] + k + sd->grid_g[0] * sd->grid_g[1]] + sd->u_now[(i * sd->grid_g[1] + j) * sd->grid_g[2] + k - sd->grid_g[0] * sd->grid_g[1]]);
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

void ApplyLaplaceFilter(Subdomain *subdomain)
{

}

void SetBoundaryConditions(Subdomain *sd)
{
    int bc = 10.0;
    if (sd->grid_l[2] == 0) // 2D problem
    {
        for (int i = 0; i < sd->grid_l[0]; i++)
        {
            sd->u_now[i * sd->grid_l[1]] = bc; // left side
            sd->u_next[i * sd->grid_l[1]] = bc;
            if (i == 0) // fill in top
            {
                for (int j = 1; j < sd->grid_l[1]; j++)
                {
                    sd->u_now[i * sd->grid_l[1] + j] = bc;
                    sd->u_next[i * sd->grid_l[1] + j] = bc;
                }
            }
            else if (i == sd->grid_l[0] - 1) // fill in bottom
            {
                for (int j = 1; j < sd->grid_l[1]; j++)
                {
                    sd->u_now[i * sd->grid_l[1] + j] = bc;
                    sd->u_next[i * sd->grid_l[1] + j] = bc;
                }
            }
            else // fill in the rest of the terms on the right side
            {
                sd->u_now[i*sd->grid_g[1] + sd->grid_g[1] - 1] = bc;
                sd->u_next[i*sd->grid_g[1] + sd->grid_g[1] - 1] = bc;
            }
        }
    }
    else // 3D problem
    { // this code doesn't work, switch to shape array
        for (int i = 0; i < sd->grid_l[0]; i++)
        {
            for (int j = 0; j < sd->grid_l[1]; j++)
            {
                for (int k = 0; k < sd->grid_l[2]; k++)
                {
                    if (i == 0 || i == sd->grid_l[0] || j == 0 || j == sd->grid_l[1] || k == 0 || k == sd->grid_l[2])
                    {
                        sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = 10;
                        sd->u_next[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = 10;
                    }
                }
            }
        }
    }

} // end void SetBoundaryConditions(Subdomain sd)

void SetInitialConditions(Subdomain *sd)
{
    int ic = 0.0;
    for (int i = 1; i < sd->grid_l[0] - 1; i++)
    {
        for (int j = 1; j < sd->grid_l[1] - 1; j++)
        {
            if (sd->grid_l[2] == 0) // 2D problem
            {
                sd->u_now[i * sd->grid_l[1] + j] = ic;
            }
            else // 3D problem
            {
                for (int k = 1; k < sd->grid_l[2] - 1; k++)
                {
                    sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k] = ic;
                }
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
    sd->dx = sd->dims_g[0] / sd->grid_g[0];
    sd->dy = sd->dims_g[1] / sd->grid_g[1];
    sd->dz = sd->dims_g[2] / sd->grid_g[2];

    // stability factors
    sd->mu_x = sd->dt / pow(sd->dx, 2.0);
    sd->mu_y = sd->dt / pow(sd->dy, 2.0);
    sd->mu_z = sd->dt / pow(sd->dz, 2.0);

    // make sure step sizes are appropriate for stable results
    if (sd->mu_x > 0.5) { fprintf(stderr, "Step sizes are too large. mu_x = %f > 0.5\n", sd->mu_x); }
    if (sd->mu_y > 0.5) { fprintf(stderr, "Step sizes are too large. mu_y = %f > 0.5\n", sd->mu_y); }
    if (sd->mu_z > 0.5) { fprintf(stderr, "Step sizes are too large. mu_z = %f > 0.5\n", sd->mu_z); }
    exit(1);
} // end void DetermineStepSizes(Subdomain sd)



