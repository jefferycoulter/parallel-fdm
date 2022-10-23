#include "parallel.h"
#include "fdm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Subdomain *CreateSubdomain(int nproc, int rank)
{   
    Subdomain *subdomain = malloc(sizeof(Subdomain));
    GetInput(subdomain, nproc, rank);
    DetermineStepSizes(subdomain);
    SplitProcessorsAlongDims(subdomain);
    DetermineSubdomainGridBounds(subdomain);
    //CreateSubdomainType(subdomain);
    //SetupCollectSubdomainData(subdomain);

    return subdomain;
} // end Subdomain *CreateSubdomain(int nproc, int rank)

void GetInput(Subdomain *sd, int n_proc, int rank)
{
    sd->n_proc = n_proc;
    sd->rank = rank;
    sd->n_dims = 2;
    sd->dims_g[0] = 40;
    sd->dims_g[1] = 40;
    sd->dims_g[2] = 0;

    sd->grid_g[0] = 200;
    sd->grid_g[1] = 200;
    sd->grid_g[2] = 1;
    sd->dt = 0.005; // 50 milliseconds
} // end void GetInput(Subdomain *sd, int n_proc, int rank)

void SplitProcessorsAlongDims(Subdomain *sd)
{
    (*sd).n_proc_dim[0] = 1;
    (*sd).n_proc_dim[1] = (*sd).n_proc;
    (*sd).n_proc_dim[2] = 1;
} // end void SplitProcessorsAlongDims(Subdomain sd)

void DetermineSubdomainGridBounds(Subdomain *sd)
{
    if (sd->n_dims == 2)
    {
        // assign the process for this subdomain cartesian coordinates.
        // to work out these equations, think of converting 10,000 seconds to hr:min:s (z:y:x)
        sd->coords[0] = floor(sd->rank / sd->n_proc_dim[1]);
        sd->coords[1] = sd->rank - sd->coords[0] * sd->n_proc_dim[1];
        sd->coords[2] = 0;

        // compute number of fdm grid cells along each dimension that will belong to a given process
        for (int n = 0; n < sd->n_dims; n++)
        {
            sd->grid_l[n] = sd->grid_g[n] / sd->n_proc_dim[n];

            // determine the bounds of the subdomain from the global perspective
            if (sd->n_proc_dim[n] == 1) // local grid is same as global grid along this dimension
            {
                // corresponds to i_start, j_start, k_start
                sd->bounds_l[2 * n] = 0;
                // corresponds to i_end, j_end, k_end
                sd->bounds_l[2 * n + 1] = sd->grid_g[n];
            }
            else // divide the global grid size amongst processes
            {
                // corresponds to i_start, j_start, k_start
                sd->bounds_l[2 * n] = sd->grid_l[n] * sd->coords[n];
                // corresponds to i_end, j_end, k_end
                sd->bounds_l[2 * n + 1] = sd->grid_l[n] * (sd->coords[n] + 1);
            }

            // this initializes the z slot to 1 in case the problem is only 2D.  if the problem is 3D
            // then these will be given the correct value on the next iteration of the for loop
            if (n == 1)
            {
                sd->grid_l[n + 1] = 1;
                sd->bounds_l[2 * n + 2] = 0;
                sd->bounds_l[2 * n + 3] = 0;
            }
        } 
    }
    else if (sd->n_dims == 3)
    {
        // assign the process for this subdomain cartesian coordinates.
        // to work out these equations, think of converting 10,000 seconds to hr:min:s (z:y:x)
        sd->coords[2] = floor(sd->rank / (sd->n_proc_dim[0] * sd->n_proc_dim[1]));
        sd->coords[1] = floor((sd->rank - sd->coords[2] * sd->n_proc_dim[0] * sd->n_proc_dim[1]) / sd->n_proc_dim[1]);
        sd->coords[0] = (sd->rank - sd->coords[2] * sd->n_proc_dim[0] * sd->n_proc_dim[1]) - sd->n_proc_dim[0] * sd->coords[1];

        // compute number of fdm grid cells along each dimension that will belong to a given process
        for (int n = 0; n < sd->n_dims; n++)
        {
            sd->grid_l[n] = sd->grid_g[n] / sd->n_proc_dim[n];
            
            // determine the bounds of the subdomain from the global perspective
            if (sd->n_proc_dim[n] == 1) // local grid is same as global grid along this dimension
            {
                // corresponds to i_start, j_start, k_start
                sd->bounds_l[2 * n] = 0;
                // corresponds to i_end, j_end, k_end
                sd->bounds_l[2 * n + 1] = sd->grid_g[n];
            }
            else // divide the global grid size amongst processes
            {
                // corresponds to i_start, j_start, k_start
                sd->bounds_l[2 * n] = sd->grid_l[n] * sd->coords[n];
                // corresponds to i_end, j_end, k_end
                sd->bounds_l[2 * n + 1] = sd->grid_l[n] * (sd->coords[n] + 1);
            }
        }
    }
    else
    {
        fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
        exit(1);
    }
} // end void DetermineSubdomainGridBounds(Subdomain *sd)

/*
void AllocateArraysFDM(Subdomain *sd)
{
    if ((*sd).n_proc == 1)
    {
        (*sd).shape_now = (int*)malloc(sizeof(int) * (*sd).grid_l[0]);
        (*sd).shape_next = (int*)malloc(sizeof(int) * (*sd).grid_l[0]);
        (*sd).u_now = (float*)malloc(sizeof(float) * (*sd).grid_l[0]);
        (*sd).u_next = (float*)malloc(sizeof(float) * (*sd).grid_l[0]);
        for (int i = 0; i < (*sd).grid_l[0]; i++)
        {
            (*sd).shape_now[i] = (int*)malloc(sizeof(int) * (*sd).grid_l[1]);
            (*sd).shape_next[i] = (int*)malloc(sizeof(int) * (*sd).grid_l[1]);
            (*sd).u_now[i] = (float*)malloc(sizeof(float) * (*sd).grid_l[1]);
            (*sd).u_next[i] = (float*)malloc(sizeof(float) * (*sd).grid_l[1]);

            for (int j = 0; j < (*sd).dims_l[1]; j++)
            {
                (*sd).shape_now[i][j] = (int*)malloc(sizeof(int) * (*sd).grid_l[2]);
                (*sd).shape_next[i][j] = (int*)malloc(sizeof(int) * (*sd).grid_l[2]);
                (*sd).u_now[i][j] = (float*)malloc(sizeof(float) * (*sd).grid_l[2]);
                (*sd).u_next[i][j] = (float*)malloc(sizeof(float) * (*sd).grid_l[2]);
            }
        }
    }
    else // add space for ghost cells
    {
        (*sd).shape_now = (int*)malloc(sizeof(int) * (*sd).grid_l[0] + 2);
        (*sd).shape_next = (int*)malloc(sizeof(int) * (*sd).grid_l[0] + 2);
        (*sd).u_now = (float*)malloc(sizeof(float) * (*sd).grid_l[0] + 2);
        (*sd).u_next = (float*)malloc(sizeof(float) * (*sd).grid_l[0] + 2);
        for (int i = 0; i < (*sd).grid_l[0]; i++)
        {
            (*sd).shape_now[i] = (int*)malloc(sizeof(int) * (*sd).grid_l[1] + 2);
            (*sd).shape_next[i] = (int*)malloc(sizeof(int) * (*sd).grid_l[1] + 2);
            (*sd).u_now[i] = (float*)malloc(sizeof(float) * (*sd).grid_l[1] + 2);
            (*sd).u_next[i] = (float*)malloc(sizeof(float) * (*sd).grid_l[1] + 2);

            for (int j = 0; j < (*sd).dims_l[1]; j++)
            {
                (*sd).shape_now[i][j] = (int*)malloc(sizeof(int) * (*sd).grid_l[2] + 2);
                (*sd).shape_next[i][j] = (int*)malloc(sizeof(int) * (*sd).grid_l[2] + 2);
                (*sd).u_now[i][j] = (float*)malloc(sizeof(float) * (*sd).grid_l[2] + 2);
                (*sd).u_next[i][j] = (float*)malloc(sizeof(float) * (*sd).grid_l[2] + 2);
            }
        }
    }

    // global shape array
    (*sd).shape_g = (int*)malloc(sizeof(int) * (*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2]);
    memset((*sd).shape_g, 0, sizeof(int) * (*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2]);

    // global solution array
    (*sd).u_g = (float*)malloc(sizeof(float) * (*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2]);
    memset((*sd).u_g, 0, sizeof(float) * (*sd).grid_g[0] * (*sd).grid_g[1] * (*sd).grid_g[2]);
} // end void AllocateArraysFDM(Subdomain sd)
 

void SubdomainCleanUp(Subdomain *sd)
{
    free((*sd).shape_now);
    free((*sd).shape_next);
    free((*sd).shape_g);
    free((*sd).u_g);
    free((*sd).u_now);
    free((*sd).u_next);
    free((*sd).send_counts);
    free((*sd).displs);

    MPI_Type_free(&(*sd).subdomain_type);
} // void SubdomainCleanUp(Subdomain *sd)
*/