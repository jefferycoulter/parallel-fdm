#include "parallel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Subdomain *CreateSubdomain(int nproc, int rank)
{   
    Subdomain *subdomain = malloc(sizeof(Subdomain));
    GetInput(subdomain, nproc, rank);
    SplitProcessorsAlongDims(subdomain);
    DetermineSubdomainGridBounds(subdomain);
    AllocateArraysFDM(subdomain);
    CreateSubdomainType(subdomain);
    SetupCollectSubdomainData(subdomain);

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

    // have to initialize this to something, otherwise mpi throws memory access errors
    sd->n_proc_dim2D[0] = sd->n_proc_dim2D[1] = 0;
    sd->n_proc_dim3D[0] = sd->n_proc_dim3D[1] = sd->n_proc_dim3D[2] = 0;

    sd->grid_g[0] = 800;
    sd->grid_g[1] = 800;
    sd->grid_g[2] = 1;
    sd->dt = 0.005; // 50 milliseconds
} // end void GetInput(Subdomain *sd, int n_proc, int rank)

void SplitProcessorsAlongDims(Subdomain *sd)
{
    // use this until i figure out how to send 3D subarrays between processors
    if (sd->n_dims == 2)
    {
        // this would divide processors evenly among dimensions
        MPI_Dims_create(sd->n_proc, sd->n_dims, sd->n_proc_dim2D);

        // split along y direction
        //sd->n_proc_dim2D[0] = sd->n_proc;
    }
    else if (sd->n_dims == 3)
    {
        // this would divide processors evenly among dimensions
        MPI_Dims_create(sd->n_proc, sd->n_dims, sd->n_proc_dim3D);
        // split along z direction
        //sd->n_proc_dim[0] = sd->n_proc;
    }
    else
    {
        fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
        exit(1);
    }
} // end void SplitProcessorsAlongDims(Subdomain sd)

void DetermineSubdomainGridBounds(Subdomain *sd)
{
    if (sd->n_dims == 2)
    {
        // assign the process for this subdomain cartesian coordinates.
        // to work out these equations, think of converting 10,000 seconds to hr:min:s (z:y:x)
        sd->coords[0] = floor(sd->rank / sd->n_proc_dim2D[1]);
        sd->coords[1] = sd->rank - sd->coords[0] * sd->n_proc_dim2D[1];
        sd->coords[2] = 0;

        // compute number of fdm grid cells along each dimension that will belong to a given process
        for (int n = 0; n < sd->n_dims; n++)
        {
            sd->grid_l[n] = sd->grid_g[n] / sd->n_proc_dim2D[n];

            // determine the bounds of the subdomain from the global perspective
            if (sd->n_proc_dim2D[n] == 1) // local grid is same as global grid along this dimension
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
        sd->coords[2] = floor(sd->rank / (sd->n_proc_dim3D[0] * sd->n_proc_dim3D[1]));
        sd->coords[1] = floor((sd->rank - sd->coords[2] * sd->n_proc_dim3D[0] * sd->n_proc_dim3D[1]) / sd->n_proc_dim3D[1]);
        sd->coords[0] = (sd->rank - sd->coords[2] * sd->n_proc_dim3D[0] * sd->n_proc_dim3D[1]) - sd->n_proc_dim3D[0] * sd->coords[1];

        // compute number of fdm grid cells along each dimension that will belong to a given process
        for (int n = 0; n < sd->n_dims; n++)
        {
            sd->grid_l[n] = sd->grid_g[n] / sd->n_proc_dim3D[n];
            
            // determine the bounds of the subdomain from the global perspective
            if (sd->n_proc_dim3D[n] == 1) // local grid is same as global grid along this dimension
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

void CreateSubdomainType(Subdomain *sd)
{
    MPI_Datatype subdomain;
    if (sd->n_dims == 2)
    {
        int sizes_g[2] = {sd->grid_g[0], sd->grid_g[1]};
        int sizes_l[2] = {sd->grid_l[0], sd->grid_l[1]};
        int starts[2] = {0, 0};

        MPI_Type_create_subarray(sd->n_dims, sizes_g, sizes_l, starts, MPI_ORDER_C, MPI_INT, &subdomain);
        MPI_Type_create_resized(subdomain, 0, sd->grid_l[1]*sizeof(int), &(*sd).subdomain_type);
        MPI_Type_commit(&(*sd).subdomain_type);
    }
    else if (sd->n_dims == 3)
    {
        int starts[3] = {0, 0, 0};
        MPI_Type_create_subarray(sd->n_dims, sd->grid_g, sd->grid_l, starts, MPI_ORDER_C, MPI_INT, &subdomain);
        MPI_Type_create_resized(subdomain, 0, sd->grid_l[0] * sizeof(int), &(*sd).subdomain_type);
        MPI_Type_commit(&(*sd).subdomain_type);
    }
    else
    {
        fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
        exit(1);
    }

    MPI_Type_free(&subdomain);
} // end void CreateSubdomainType(Subdomain *sd)

void CreateSubdomainType2(Subdomain *sd)
{
    MPI_Datatype one_slice, two_slice;
    if (sd->n_dims == 2)
    {
        MPI_Type_vector((*sd).grid_l[0], 1, sizeof(int), MPI_INT, &one_slice);
        MPI_Type_hvector((*sd).grid_l[1], 1, (*sd).grid_g[0] * sizeof(int), one_slice, &(*sd).subdomain_type);
        MPI_Type_commit(&(*sd).subdomain_type);

        MPI_Type_free(&one_slice);
    }
    else if (sd->n_dims == 3)
    {
        MPI_Type_vector((*sd).grid_l[0], 1, sizeof(int), MPI_INT, &one_slice);
        MPI_Type_hvector((*sd).grid_l[1], 1, (*sd).grid_g[0] * sizeof(int), one_slice, &two_slice);
        MPI_Type_hvector((*sd).grid_l[2], 1, (*sd).grid_g[0] * (*sd).grid_g[1] * sizeof(int), two_slice, &(*sd).subdomain_type);
        MPI_Type_commit(&(*sd).subdomain_type);

        MPI_Type_free(&one_slice);
        MPI_Type_free(&two_slice);
    }
    else
    {
        fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
        exit(1);
    }
} // void CreateSubdomainType2(Subdomain *sd)

void CollectSubdomainData(Subdomain *sd)
{
    // should use this call once 3D subdomain divisions are figured out
    MPI_Gatherv((*sd).shape_now, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, 
                (*sd).shape_g, (*sd).send_counts, (*sd).displs, (*sd).subdomain_type,
                ROOT, MPI_COMM_WORLD);
               
    //MPI_Gather((*sd).shape_next, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, 
      //         (*sd).shape_g, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT,
        //       ROOT, MPI_COMM_WORLD);

    //MPI_Gather((*sd).shape_now, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, 
      //         (*sd).shape_g, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT,
        //       ROOT, MPI_COMM_WORLD);
} // end void CollectSubdomainData(Subdomain *subdomain)

void SetupCollectSubdomainData(Subdomain *sd)
{
    if (sd->n_dims == 2)
    {
        sd->send_counts = (int*)malloc((*sd).n_proc_dim2D[0] * (*sd).n_proc_dim2D[1] * sizeof(int));
        sd->displs = (int*)malloc((*sd).n_proc_dim2D[0] * (*sd).n_proc_dim2D[1] * sizeof(int));
        memset((*sd).send_counts, 1, (*sd).n_proc_dim2D[0] * (*sd).n_proc_dim2D[1] * sizeof(int));

        int disp = 0;
        for (int i = 0; i < sd->n_proc_dim2D[0]; i++)
        {
            for (int j = 0; j < sd->n_proc_dim2D[1]; j++)
            {
                sd->displs[i * sd->n_proc_dim2D[1] + j] = disp;
                disp += 1;
            }
            disp += (sd->grid_l[0] - 1) * sd->n_proc_dim2D[1]; // 1,1 on np2 works
        }
    }
    else if (sd->n_dims == 3)
    {
        sd->send_counts = (int*)malloc((*sd).n_proc_dim3D[0] * (*sd).n_proc_dim3D[1] * (*sd).n_proc_dim3D[2] * sizeof(int));
        sd->displs = (int*)malloc((*sd).n_proc_dim3D[0] * (*sd).n_proc_dim3D[1] * (*sd).n_proc_dim3D[2] * sizeof(int));
        memset((*sd).send_counts, 1, (*sd).n_proc_dim3D[0] * (*sd).n_proc_dim3D[1] * (*sd).n_proc_dim3D[2] * sizeof(int));

        int disp = 0;
        for (int i = 0; i < sd->n_proc_dim3D[0]; i++)
        {
            for (int j = 0; j < sd->n_proc_dim3D[1]; j++)
            {
                for (int k = 0; k < sd->n_proc_dim3D[2]; k++)
                {
                    sd->displs[(i * sd->n_proc_dim3D[1] + j) * sd->n_proc_dim3D[2] + k] = disp;
                    disp += 1;
                }
                disp += (sd->grid_l[1] - 1) * sd->n_proc_dim3D[1];
            }
            disp += (sd->grid_l[1] - 1) * sd->n_proc_dim3D[1];
        }
    }
    else
    {
        fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
        exit(1);
    }
} // end void SetupCollectSubdomainData(Subdomain *subdomain)

void AllocateArraysFDM(Subdomain *sd)
{
    // size of memory to allocate
    int local_size, global_size;

    if (sd->n_dims == 2)
    {
        if (sd->n_proc == 1) // sequential, so don't need to add space for ghost cells
        {
            local_size = ((*sd).grid_l[0]) * ((*sd).grid_l[1]);
            global_size = ((*sd).grid_g[0]) * ((*sd).grid_g[1]);

            // local shape arrays
            (*sd).shape_now = (int*)malloc(sizeof(int) * local_size);
            (*sd).shape_next = (int*)malloc(sizeof(int) * local_size);
            memset((*sd).shape_now, 0, sizeof(int) * local_size);
            memset((*sd).shape_next, 0, sizeof(int) * local_size);
            
            // global shape array
            (*sd).shape_g = (int*)malloc(sizeof(int) * global_size);
            memset((*sd).shape_g, 0, sizeof(int) * global_size);

            // local solution arrays
            (*sd).u_next = (float*)malloc(sizeof(float) * local_size);
            (*sd).u_now = (float*)malloc(sizeof(float) * local_size);
            memset((*sd).u_next, 0, sizeof(float) * local_size);
            memset((*sd).u_now, 0, sizeof(float) * local_size);

            // global solution array
            (*sd).u_global = (float*)malloc(sizeof(float) * global_size);
            memset((*sd).u_global, 0, sizeof(float) * global_size);
        }
        else // parallel, so need to add space for ghost cells
        {
            local_size = ((*sd).grid_l[0] + 2) * ((*sd).grid_l[1] + 2);
            global_size = ((*sd).grid_g[0]) * ((*sd).grid_g[1]);

            // local shape arrays
            (*sd).shape_now = (int*)malloc(sizeof(int) * local_size);
            (*sd).shape_next = (int*)malloc(sizeof(int) * local_size);
            memset((*sd).shape_now, 0, sizeof(int) * local_size);
            memset((*sd).shape_next, 0, sizeof(int) * local_size);
            
            // global shape array
            (*sd).shape_g = (int*)malloc(sizeof(int) * global_size);
            memset((*sd).shape_g, 0, sizeof(int) * global_size);

            // local solution arrays
            (*sd).u_next = (float*)malloc(sizeof(float) * local_size);
            (*sd).u_now = (float*)malloc(sizeof(float) * local_size);
            memset((*sd).u_next, 0, sizeof(float) * local_size);
            memset((*sd).u_now, 0, sizeof(float) * local_size);

            // global solution array
            (*sd).u_global = (float*)malloc(sizeof(float) * global_size);
            memset((*sd).u_global, 0, sizeof(float) * global_size);
        } // end if-else for n_proc
    }
    else if (sd->n_dims == 3)
    {
        if (sd->n_proc == 1)
        {
            local_size = ((*sd).grid_l[0]) * ((*sd).grid_l[1]) * ((*sd).grid_l[2]);
            global_size = ((*sd).grid_g[0]) * ((*sd).grid_g[1]) * ((*sd).grid_g[2]);

            // local shape arrays
            (*sd).shape_now = (int*)malloc(sizeof(int) * local_size);
            (*sd).shape_next = (int*)malloc(sizeof(int) * local_size);
            memset((*sd).shape_now, 0, sizeof(int) * local_size);
            memset((*sd).shape_next, 0, sizeof(int) * local_size);
            
            // global shape array
            (*sd).shape_g = (int*)malloc(sizeof(int) * global_size);
            memset((*sd).shape_g, 0, sizeof(int) * global_size);

            // local solution arrays
            (*sd).u_next = (float*)malloc(sizeof(float) * local_size);
            (*sd).u_now = (float*)malloc(sizeof(float) * local_size);
            memset((*sd).u_next, 0, sizeof(float) * local_size);
            memset((*sd).u_now, 0, sizeof(float) * local_size);

            // global solution array
            (*sd).u_global = (float*)malloc(sizeof(float) * global_size);
            memset((*sd).u_global, 0, sizeof(float) * global_size);
        }
        else
        {
            local_size = ((*sd).grid_l[0] + 2) * ((*sd).grid_l[1] + 2) * ((*sd).grid_l[2] + 2);
            global_size = ((*sd).grid_g[0]) * ((*sd).grid_g[1]) * ((*sd).grid_g[2]);

            // local shape arrays
            (*sd).shape_now = (int*)malloc(sizeof(int) * local_size);
            (*sd).shape_next = (int*)malloc(sizeof(int) * local_size);
            memset((*sd).shape_now, 0, sizeof(int) * local_size);
            memset((*sd).shape_next, 0, sizeof(int) * local_size);
            
            // global shape array
            (*sd).shape_g = (int*)malloc(sizeof(int) * global_size);
            memset((*sd).shape_g, 0, sizeof(int) * global_size);

            // local solution arrays
            (*sd).u_next = (float*)malloc(sizeof(float) * local_size);
            (*sd).u_now = (float*)malloc(sizeof(float) * local_size);
            memset((*sd).u_next, 0, sizeof(float) * local_size);
            memset((*sd).u_now, 0, sizeof(float) * local_size);

            // global solution array
            (*sd).u_global = (float*)malloc(sizeof(float) * global_size);
            memset((*sd).u_global, 0, sizeof(float) * global_size);
        } // end if-else for n_proc
    }
    else
    {
        fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
        exit(1);
    } // end if-else for n_dims
} // end void AllocateArraysFDM(Subdomain sd)

void SubdomainCleanUp(Subdomain *sd)
{
    free((*sd).shape_now);
    free((*sd).shape_next);
    free((*sd).shape_g);
    free((*sd).u_global);
    free((*sd).u_now);
    free((*sd).u_next);
    free((*sd).send_counts);
    free((*sd).displs);

    MPI_Type_free(&(*sd).subdomain_type);
} // void SubdomainCleanUp(Subdomain *sd)

void CoordShift(Subdomain *sd, float r)
{
    //fprintf(stdout, "coord shift\n");
    int local_id; // the computations below are shifted to global coordinates, this converts back to local index
    float r_shifted; // squared length of shifted position vector of a fdm grid_g cell
    for (int i = sd->bounds_l[0]; i < sd->bounds_l[1]; i ++)
    {
        for (int j = sd->bounds_l[2]; j < sd->bounds_l[3]; j++)
        {
            if (sd->n_dims == 2) // 2D problem
            {
                local_id = (i - sd->bounds_l[0]) * sd->grid_l[1] + j - sd->bounds_l[2];
                r_shifted = pow(i - (sd->grid_g[0] / 2), 2.0) + pow(j - (sd->grid_g[1] / 2), 2.0);
                if (r_shifted <= r*r) // inside of circle
                {
                    fprintf(stdout, "rank %d: r_shifted = %f\n", sd->rank, r_shifted);
                    fprintf(stdout, "rank %d: sd->n_proc_dim[0] = %d\n", sd->rank, sd->n_proc_dim2D[0]);
                    fprintf(stdout, "rank %d: sd->n_proc_dim[1] = %d\n", sd->rank, sd->n_proc_dim2D[1]);
                    fprintf(stdout, "coordinates (%d, %d)\n", i, j);
                    fprintf(stdout, "rank %d: sd->coords[0] = %d\n", sd->rank, sd->coords[0]);
                    fprintf(stdout, "rank %d: sd->coords[1] = %d\n", sd->rank, sd->coords[1]);
                    fprintf(stdout, "rank %d: x bounds [%d, %d]\n", sd->rank, sd->bounds_l[0], sd->bounds_l[1]);
                    fprintf(stdout, "rank %d: y bounds [%d, %d]\n", sd->rank, sd->bounds_l[2], sd->bounds_l[3]);
                    sd->shape_now[local_id] = 1;
                }
                else // outside of circle
                {
                    sd->shape_now[local_id] = 0;
                }
            }
            else if (sd->n_dims == 3)// 3D problem
            {
                for (int k = sd->bounds_l[4]; k < sd->bounds_l[5]; k++)
                {
                    local_id = ((i - sd->bounds_l[0]) * sd->grid_l[1] + j - sd->bounds_l[2]) * sd->grid_l[2] + k - sd->bounds_l[4];
                    r_shifted = pow(i - (sd->grid_g[0] / 2), 2.0) + \
                                pow(j - (sd->grid_g[1] / 2), 2.0) + \
                                pow(k - (sd->grid_g[2] / 2), 2.0);
                    if (r_shifted <= r*r) // inside of circle
                    {
                        sd->shape_now[local_id] = 1;
                    }
                    else // outside of circle
                    {
                        sd->shape_now[local_id] = 0;
                    }
                } // end k loop
            } // end if-else
            else
            {
                fprintf(stderr, "Invalid number of dimensions: %d.  Need to choose either 2D or 3D.\n", sd->n_dims);
                exit(1);
            }
        } // end j loop
    } // end i loop

    MPI_Barrier(MPI_COMM_WORLD);
} // end void CoordShift(Subdomain sd, float r)