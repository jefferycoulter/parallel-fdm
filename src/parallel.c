#include "parallel.h"

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Subdomain *CreateSubdomain(int nproc, int rank)
{   
    Subdomain *subdomain = malloc(sizeof(Subdomain));
    GetInput(subdomain, nproc, rank);
    fprintf(stdout, "sd->rank = %d\n", subdomain->rank);
    fprintf(stdout, "sd->n_dims = %d\n", subdomain->n_dims);
    SplitProcessorsAlongDims(subdomain);
    fprintf(stdout, "sd->n_proc_dim[0] = %d\n", subdomain->n_proc_dim[0]);
    fprintf(stdout, "sd->n_proc_dim[1] = %d\n", subdomain->n_proc_dim[1]);
    fprintf(stdout, "sd->n_proc_dim[2] = %d\n", subdomain->n_proc_dim[2]);
    GetSubdomainGridBounds(subdomain);
    ComputeSubdomainGrid(subdomain);
    AllocateArraysFDM(subdomain);

    return subdomain;
}

void GetInput(Subdomain *sd, int n_proc, int rank)
{
    sd->n_proc = n_proc;
    sd->rank = rank,
    sd->n_dims = 2;
    sd->dims_g[0] = 20;
    sd->dims_g[1] = 20;
    sd->dims_g[2] = 0;

    sd->n_proc_dim[0] = 0;
    sd->n_proc_dim[1] = 0;
    sd->n_proc_dim[2] = 0;

    sd->grid_g[0] = 200;
    sd->grid_g[1] = 200;
    sd->grid_g[2] = 1;
    sd->dt = 0.05; // 50 milliseconds
}

void SplitProcessorsAlongDims(Subdomain *sd)
{
    MPI_Dims_create(sd->n_proc, sd->n_dims, sd->n_proc_dim);
} // end void SplitProcessorsAlongDims(Subdomain sd)

void GetSubdomainGridBounds(Subdomain *sd)
{
    for (int n = 0; n < (*sd).n_dims; n++)
    {
        if ((*sd).n_proc_dim[n] == 1) // local grid is same as global grid along this dimension
        {
            // corresponds to i_start, j_start, k_start
            (*sd).bounds_l[2 * n] = 0;
            // corresponds to i_end, j_end, k_end
            (*sd).bounds_l[2 * n + 1] = (*sd).grid_g[n];
        }
        else // divide the global grid size amongst processes
        {
            // corresponds to i_start, j_start, k_start
            (*sd).bounds_l[2 * n] = (*sd).rank * ((*sd).grid_g[n] / (*sd).n_proc_dim[n]);
            // corresponds to i_end, j_end, k_end
            (*sd).bounds_l[2 * n + 1] = (*sd).bounds_l[2 * n] + ((*sd).grid_g[n] / (*sd).n_proc_dim[n]);
        }

        // this zeros out the z slots in case the problem is only 2D.  if the problem is 3D
        // then these will be given the correct value on the next iteration of the for loop
        if (n == 1)
        {
            (*sd).bounds_l[2 * n + 2] = 0;
            (*sd).bounds_l[2 * n + 3] = 0;
        }
    }
} // end void GetSubdomainGridBounds(Subdomain sd)

void ComputeSubdomainGrid(Subdomain *sd)
{
    // loop over dimensions
    for(int n = 0; n < (*sd).n_dims; n++)
    {

        (*sd).grid_l[n] = (*sd).bounds_l[2 * n + 1] - (*sd).bounds_l[2 * n];

        // this zeros out the z slot in case the problem is only 2D.  if the problem is 3D
        // then this will be given the correct value on the next iteration of the for loop
        if (n == 1)
        {
            (*sd).grid_l[n + 1] = 1;
        }
    }
} // end void ComputeSubdomainDims(Subdomain sd)

void AllocateArraysFDM(Subdomain *sd)
{
    // note:
    // sd.grid_g[2] = sd.grid_l[2] = 1 for 2D case, so multiplying it here doesn't matter 
    // if the problem is 2D

    // shape arrays
    (*sd).shape_arr_l = (int*)malloc((*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof(int));
    (*sd).shape_arr_g = (int*)malloc((*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2] * sizeof(int));
    memset((*sd).shape_arr_l, 0, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof(int));
    memset((*sd).shape_arr_g, 0, (*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2] * sizeof(int));

    // local solution arrays
    (*sd).u_next = (float*)malloc((*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof(float));
    (*sd).u_now = (float*)malloc((*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof(float));
    memset((*sd).u_next, 0, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof(((*sd).u_next)));
    memset((*sd).u_now, 0, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2] * sizeof(((*sd).u_now)));

    // global solution array
    (*sd).u_global = (float*)malloc((*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2] * sizeof(float));
    memset((*sd).u_global, 0, (*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2] * sizeof(float));
} // end void AllocateArraysFDM(Subdomain sd)

void CoordShift(Subdomain *sd, float r)
{
    float r_shifted; // squared length of shifted position vector of a fdm grid_g cell
    for (int i = (*sd).bounds_l[0]; i < (*sd).bounds_l[1]; i ++)
    {
        for (int j = (*sd).bounds_l[2]; j < (*sd).bounds_l[3]; j++)
        {
            if ((*sd).n_dims == 2) // 2D problem
            {
                r_shifted = pow(i - ((*sd).grid_g[0] / 2), 2.0) + pow(j - ((*sd).grid_g[1] / 2), 2.0);
                if ((int)r_shifted <= (int)r*r) // inside of circle
                {
                    (*sd).shape_arr_l[(i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]] = 1;
                }
                else // outside of circle
                {
                    (*sd).shape_arr_l[(i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]] = 0;
                }
            }
            else // 3D problem
            {
                for (int k = (*sd).bounds_l[4]; k < (*sd).bounds_l[5]; k++)
                {
                    r_shifted = pow(i - ((*sd).grid_g[0] / 2), 2.0) + \
                                pow(j - ((*sd).grid_g[1] / 2), 2.0) + \
                                pow(k - ((*sd).grid_g[2] / 2), 2.0);

                    if ((int)r_shifted <= (int)r*r) // inside of circle
                    {
                        (*sd).shape_arr_l[((i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]) * (*sd).grid_l[2] + k - (*sd).bounds_l[4]] = 1;
                    }
                    else // outside of circle
                    {
                        (*sd).shape_arr_l[((i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]) * (*sd).grid_l[2] + k - (*sd).bounds_l[4]] = 0;
                    }
                } // end k loop
            } // end if-else
        } // end j loop
    } // end i loop
} // end void CoordShift(Subdomain sd, float r)