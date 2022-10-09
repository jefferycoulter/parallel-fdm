#include "parallel.h"

#include <mpi.h>
#include <math.h>
#include <stdio.h>

void PrepareSubdomains(Subdomain *subdomain)
{
    SplitProcessorsAlongDims(subdomain);
    GetSubdomainGridBounds(subdomain);
    ComputeSubdomainGrid(subdomain);
} // end void PrepareSubdomains(Subdomain subdomain)

void SplitProcessorsAlongDims(Subdomain *sd)
{
    MPI_Dims_create((*sd).n_proc, (*sd).n_dims, (*sd).n_proc_dim);
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

        fprintf(stdout, "rank %d: global grid along dim %d is %d\n", (*sd).rank, n, (*sd).grid_g[n]);
        fprintf(stdout, "rank %d: n_proc_dim[%d] is %d\n", (*sd).rank, n, (*sd).n_proc_dim[n]);
        fprintf(stdout, "rank %d: lower grid id along dim %d is %d, upper is %d\n", (*sd).rank, n, (*sd).bounds_l[2 * n], (*sd).bounds_l[2 * n + 1]);
        // this zeros out the z slots in case the problem is only 2D.  if the problem is 3D
        // then these will be given the correct value on the next iteration of the for loop
        if (n == 2)
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
        fprintf(stdout, "rank %d: grid_l[%d] is %d\n", (*sd).rank, n, (*sd).bounds_l[2 * n + 1] - (*sd).bounds_l[2 * n]);

        // this zeros out the z slot in case the problem is only 2D.  if the problem is 3D
        // then this will be given the correct value on the next iteration of the for loop
        if (n == 1)
        {
            (*sd).grid_l[n + 1] = 1;
        }
    }
} // end void ComputeSubdomainDims(Subdomain sd)

void CoordShift(Subdomain *sd, float r)
{

    //fprintf(stdout, "lower grid id along dim 0 is %d, upper is %d\n", sd.bounds_l[0], sd.bounds_l[1]);
    //fprintf(stdout, "lower grid id along dim 1 is %d, upper is %d\n", sd.bounds_l[2], sd.bounds_l[3]);
    float r_shifted; // squared length of shifted position vector of a fdm grid_g cell
    for (int i = (*sd).bounds_l[0]; i < (*sd).bounds_l[1]-1; i ++)
    {
        for (int j = (*sd).bounds_l[2]; j < (*sd).bounds_l[3]-1; j++)
        {
            if ((*sd).n_dims == 2) // 2D problem
            {
                r_shifted = pow(i - ((*sd).grid_g[0] / 2), 2.0) + pow(j - ((*sd).grid_g[1] / 2), 2.0);
                if ((int)r_shifted <= (int)r*r) // inside of circle
                {
                    fprintf(stdout, "rank %d: index is %d\n", (*sd).rank, (i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]);
                    fprintf(stdout, "rank %d:max index is %d\n", (*sd).rank, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2]);
                    (*sd).shape_arr_l[(i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]] = 1;
                    //fprintf(stdout, "updated shape array\n");
                }
                else // outside of circle
                {
                    fprintf(stdout, "rank %d: index is %d\n", (*sd).rank, (i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]);
                    fprintf(stdout, "rank %d:max index is %d\n", (*sd).rank, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2]);
                    (*sd).shape_arr_l[(i - (*sd).bounds_l[0]) * (*sd).grid_l[1] + j - (*sd).bounds_l[2]] = 0;
                    //fprintf(stdout, "didn't update shape array\n");
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