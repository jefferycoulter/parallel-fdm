#include "fdm.h"
#include "io.h"

#include <mpi.h>

#include <stdlib.h>
#include <math.h>
#include <memory.h>

void ComputeFD(Subdomain sd, int bc, int time_steps)
{
    WriteData(sd.u_now, sd.grid); // save initial conditions
    for (int t = 1; t < time_steps+1; t ++)
    {
        FTCS(sd, bc); // compute difference
       // WriteData(data.next, cells); // save each iteration
    }
}

void FTCS(Subdomain sd, int bc)
{
    InteriorFD(sd);
    BoundaryFD(sd, bc);
    // compute finite difference
    switch (bc)
    {
        case Dirichlet:
        {
            break;
        } // end case Dirichlet
        case VonNeumann:
        {
            float bc = 0.0;
            for (int x = 1; x < sd.grid[0]; x++)
            {
                for (int y = 1; y < sd.grid[1]; y++)
                {
                    if (sd.grid[2] == 0.0) // 2D problem
                    {
                        sd.u_next[x * sd.grid[1] + y] = (1 - 2 * sd.mu_x - 2 * sd.mu_y) * sd.u_now[x * sd.grid[1] + y]          \
                            + sd.mu_x * (sd.u_now[x * sd.grid[1] + y + sd.grid[0]] + sd.u_now[x * sd.grid[1] + y - sd.grid[0]])  \
                            + sd.mu_y * (sd.u_now[x * sd.grid[1] + y + 1] + sd.u_now[x * sd.grid[1] + y - 1]);
                    }
                    else // 3d problem
                    {
                        for (int z = 0; z < sd.grid[2]; z++)
                        {
                            sd.u_next[(x * sd.grid[1] + y) * sd.grid[2] + z] = (1 - 2 * sd.mu_x - 2 * sd.mu_y - 2 * sd.mu_z) * sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z]                   \
                            + sd.mu_x * (sd.u_now[(x * sd.grid[1] + y + sd.grid[0]) * sd.grid[2] + z] + sd.u_now[(x * sd.grid[1] + y - sd.grid[0]) * sd.grid[2] + z])    \
                            + sd.mu_y * (sd.u_now[(x * sd.grid[1] + y + 1) * sd.grid[2] + z] + sd.u_now[(x * sd.grid[1] + y - 1) * sd.grid[2] + z])                  \
                            + sd.mu_z * (sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z + sd.grid[0] * sd.grid[1]] + sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z - sd.grid[0] * sd.grid[1]]);
                        } // end z loop
                    } // end if-else
                } // end y loop
            } // end x loop
            break;
        } // end case VonNeumann
    } // end switch

    // copy new data to old data array for next iteration
    if (sd.n_dims == 2)
    {
        memcpy(sd.u_now, sd.u_next, sd.grid[0] * sd.grid[1] * sizeof((*sd.u_now)));
    }
    else
    {
        memcpy(sd.u_now, sd.u_next, sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof((*sd.u_now)));
    }
}

void InteriorFD(Subdomain sd)
{
    for (int x = 1; x < sd.grid[0] - 1; x++)
    {
        for (int y = 1; y < sd.grid[1] - 1; y++)
        {
            if (sd.grid[2] == 0.0) // 2D problem
            {
                sd.u_next[x * sd.grid[1] + y] = (1 - 2 * sd.mu_x - 2 * sd.mu_y) * sd.u_now[x * sd.grid[1] + y]          \
                    + sd.mu_x * (sd.u_now[x * sd.grid[1] + y + sd.grid[0]] + sd.u_now[x * sd.grid[1] + y - sd.grid[0]])  \
                    + sd.mu_y * (sd.u_now[x * sd.grid[1] + y + 1] + sd.u_now[x * sd.grid[1] + y - 1]);
            }
            else
            {
                for (int z = 1; z < sd.grid[2] - 1; z++)
                {
                    sd.u_next[(x * sd.grid[1] + y) * sd.grid[2] + z] = (1 - 2 * sd.mu_x - 2 * sd.mu_y - 2 * sd.mu_z) * sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z]                   \
                    + sd.mu_x * (sd.u_now[(x * sd.grid[1] + y + sd.grid[0]) * sd.grid[2] + z] + sd.u_now[(x * sd.grid[1] + y - sd.grid[0]) * sd.grid[2] + z])    \
                    + sd.mu_y * (sd.u_now[(x * sd.grid[1] + y + 1) * sd.grid[2] + z] + sd.u_now[(x * sd.grid[1] + y - 1) * sd.grid[2] + z])                  \
                    + sd.mu_z * (sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z + sd.grid[0] * sd.grid[1]] + sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z - sd.grid[0] * sd.grid[1]]);
                } // end z loop
            } // end if-else
        } // end y loop
    } // end x loop
}

void BoundaryFD(Subdomain sd, int bc)
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
            for (int x = 0; x < sd.grid[0]; x++)
            {
                for (int y = 0; y < sd.grid[1]; y++)
                {
                    if (sd.grid[2] == 0)
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
                        for (int z = 0; z < sd.grid[2]; z++)
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
                        }
                    }
                    break;
                }
            }
            
        }
    }
}

void SetBoundaryConditions(Subdomain sd)
{
    int bc = 10.0;
    if (sd.grid[2] == 0.0) // 2D problem
    {
        for (int x = 0; x < sd.grid[0]; x++)
        {
            sd.u_now[x * sd.grid[1]] = bc; // left side
            sd.u_next[x * sd.grid[1]] = bc;
            if (x == 0) // fill in top
            {
                for (int y = 1; y < sd.grid[1]; y++)
                {
                    sd.u_now[x * sd.grid[1] + y] = bc;
                    sd.u_next[x * sd.grid[1] + y] = bc;
                }
            }
            else if (x == sd.grid[0] - 1) // fill in bottom
            {
                for (int y = 1; y < sd.grid[1]; y++)
                {
                    sd.u_now[x * sd.grid[1] + y] = bc;
                    sd.u_next[x * sd.grid[1] + y] = bc;
                }
            }
            else // fill in the rest of the terms on the right side
            {
                sd.u_now[x*sd.grid[1] + sd.grid[1] - 1] = bc;
                sd.u_next[x*sd.grid[1] + sd.grid[1] - 1] = bc;
            }
        }
    }
    else // 3D problem
    {
        for (int i = 0; i < sd.grid[0]; i++)
        {
            for (int j = 0; j < sd.grid[1]; j++)
            {
                for (int k = 0; k < sd.grid[2]; k++)
                {
                    if (i == 0 || i == sd.grid[0] || j == 0 || j == sd.grid[1] || k == 0 || k == sd.grid[2])
                    {
                        sd.u_now[(i * sd.grid[1] + j) * sd.grid[2] + k] = 10;
                        sd.u_next[(i * sd.grid[1] + j) * sd.grid[2] + k] = 10;
                    }
                }
            }
        }
    }

}

void SetInitialConditions(Subdomain sd)
{
    int ic = 0.0;
    for (int x = 1; x < sd.grid[0] - 1; x++)
    {
        for (int y = 1; y < sd.grid[1] - 1; y++)
        {
            if (sd.grid[2] == 0) // 2D problem
            {
                sd.u_now[x * sd.grid[1] + y] = ic;
            }
            else // 3D problem
            {
                for (int z = 1; z < sd.grid[2] - 1; z++)
                {
                    sd.u_now[(x * sd.grid[1] + y) * sd.grid[2] + z] = ic;
                }
            }
        }
    }
}

void DiscretizeSubdomain(Subdomain sd)
{
    if(sd.grid[2] == 0)
    {
        sd.grid[0] = sd.dims_l[0] / sd.dx;
        sd.grid[1] = sd.dims_l[1] / sd.dy;
        sd.grid[2] = 0;
    }
    else
    {
        sd.grid[0] = sd.dims_l[0] / sd.dx;
        sd.grid[1] = sd.dims_l[1] / sd.dy;
        sd.grid[2] = sd.dims_l[2] / sd.dz;
    }
}

void GenerateArraysFDM(Subdomain sd)
{
    if (sd.grid[2] == 0.0) // 2D problem
    {
        sd.shape_arr = (float*)malloc(sd.grid[0] * sd.grid[1] * sizeof(float));
        sd.u_next = (float*)malloc(sd.grid[0] * sd.grid[1] * sizeof(float));
        sd.u_now = (float*)malloc(sd.grid[0] * sd.grid[1] * sizeof(float));
        memset(sd.shape_arr, 0, sd.grid[0] * sd.grid[1] * sizeof((*sd.u_next)));
        memset(sd.u_next, 0, sd.grid[0] * sd.grid[1] * sizeof((*sd.u_next)));
        memset(sd.u_now, 0, sd.grid[0] * sd.grid[1] * sizeof((*sd.u_now)));
    }
    else // 3D problem
    {
        sd.shape_arr = (float*)malloc(sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof(float));
        sd.u_next = (float*)malloc(sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof(float));
        sd.u_now = (float*)malloc(sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof(float));
        memset(sd.shape_arr, 0, sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof((*sd.u_next)));
        memset(sd.u_next, 0, sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof((*sd.u_next)));
        memset(sd.u_now, 0, sd.grid[0] * sd.grid[1] * sd.grid[2] * sizeof((*sd.u_now)));
    }
}



