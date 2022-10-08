#include "parallel.h"

#include <math.h>

void GetSubdomainBounds(Subdomain sd)
{
    // need to probably add something to adjust global size if dims_g / n isn't an integer
    for (int n = 0; n < sd.n_dims; n++)
    {
        // corresponds to x_start, y_start, z_start
        sd.bounds_l[2 * n] = (int)floor(sd.rank * (sd.dims_g[n] / sd.n_proc));
        // corresponds to x_end, y_end, z_end
        sd.bounds_l[2 * n + 1] = (int)floor(sd.bounds_l[n] + (sd.dims_g[n] / sd.n_proc));

        // this zeros out the z slots in case the problem is only 2D.  if the problem is 3D
        // then these will be given the correct value on the next iteration of the for loop
        if (n == 2)
        {
            sd.bounds_l[2 * n + 2] = 0;
            sd.bounds_l[2 * n + 3] = 0;
        }
    }
}

void ComputeSubdomainDims(Subdomain sd)
{
    for(int n = 0; n < sd.n_dims; n++)
    {
        sd.dims_l[n] = (int)floor(sd.dims_g[n] / sd.n_proc);
        // this zeros out the z slot in case the problem is only 2D.  if the problem is 3D
        // then this will be given the correct value on the next iteration of the for loop
        if (n == 2)
        {
            sd.dims_l[n + 1] = 0;
        }
    }
}

void CoordShift(Subdomain sd, float r)
{
    for (int i = sd.bounds_l[0]; i < sd.bounds_l[1]; i ++)
    {
        for (int j = sd.bounds_l[2]; j < sd.bounds_l[3]; j++)
        {
            if (sd.n_dims ==2) // 2D problem
            {
                float r_shift = pow(i - (sd.dims_g[0] / 2), 2.0) + pow(j - (sd.dims_g[1] / 2), 2.0);
                if ((int)r_shift <= (int)r*r) // inside of circle
                {
                    sd.shape_arr[i * sd.grid[1] + j] = 1;
                }
                else // outside of circle
                {
                    sd.shape_arr[i * sd.grid[1] + j] = 0;
                }
            }
            else // 3D problem
            {
                for (int k = sd.bounds_l[4]; k < sd.bounds_l[5]; k++)
                {
                    float r_shift = pow(i - (sd.dims_g[0] / 2), 2.0) + \
                                    pow(j - (sd.dims_g[1] / 2), 2.0) + \
                                    pow(k - (sd.dims_g[2] / 2), 2.0);

                    if ((int)r_shift <= (int)r*r) // inside of circle
                    {
                        sd.shape_arr[(i * sd.grid[1] + j) * sd.grid[2] + k] = 1;
                    }
                    else // outside of circle
                    {
                        sd.shape_arr[(i * sd.grid[1] + j) * sd.grid[2] + k] = 0;
                    }
                }
            }
        }
    }
}