#include "parallel.h"
#include "io.h"
#include "fdm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Subdomain *CreateSubdomain(int nproc, int rank)
{   
    Subdomain *subdomain = malloc(sizeof(Subdomain));
    GetInput(subdomain, nproc, rank);
    SplitProcessorsAlongDims(subdomain);
    GetNeighbors(subdomain);
    DetermineSubdomainGridBounds(subdomain);
    AllocateArraysFDM(subdomain);

    return subdomain;
} // end Subdomain *CreateSubdomain(int nproc, int rank)

void GetInput(Subdomain *sd, int n_proc, int rank)
{
    sd->n_proc = n_proc;
    sd->rank = rank;
    sd->n_dims = 3;
    sd->dims_g[0] = 40;
    sd->dims_g[1] = 40;
    sd->dims_g[2] = 40;

    sd->grid_g[0] = 100;
    sd->grid_g[1] = 100;
    sd->grid_g[2] = 100;
    sd->dt = 0.005; // 50 milliseconds
} // end void GetInput(Subdomain *sd, int n_proc, int rank)

void SplitProcessorsAlongDims(Subdomain *sd)
{
    // use this until i figure out how to send 3D subarrays between processors
    if (sd->n_dims == 2)
    {
        // split along y direction
        sd->n_proc_dim[0] = sd->n_proc;
        sd->n_proc_dim[1] = 1;
    }
    else if (sd->n_dims == 3)
    {
        // split along y direction
        sd->n_proc_dim[0] = sd->n_proc;
        sd->n_proc_dim[1] = 1;
        sd->n_proc_dim[2] = 1;
    }
} // end void SplitProcessorsAlongDims(Subdomain sd)

void GetNeighbors(Subdomain *sd)
{
    sd->neighbors[0] = MPI_PROC_NULL;
    sd->neighbors[1] = MPI_PROC_NULL;
    
    if (sd->n_dims ==2 )
    { 
        int periods[2] = {0, 0}; 
        int proc_dim[2] = {sd->n_proc_dim[0], sd->n_proc_dim[1]};
        MPI_Cart_create(MPI_COMM_WORLD, sd->n_dims, proc_dim, periods, 0, &(*sd).COMM_FDM);
        MPI_Cart_shift(sd->COMM_FDM, 0, 1, &(*sd).neighbors[Down], &(*sd).neighbors[Up]);
    }
    else if (sd->n_dims == 3)
    {
        int periods[3] = {0, 0, 0}; 
        MPI_Cart_create(MPI_COMM_WORLD, sd->n_dims, sd->n_proc_dim, periods, 0, &(*sd).COMM_FDM);
        MPI_Cart_shift(sd->COMM_FDM, 0, 1, &(*sd).neighbors[Down], &(*sd).neighbors[Up]);
    }
} // end void GetNeighbors(Subdomain *sd)

void ShareGhosts(Subdomain *sd, int type)
{

    int id_send_up = (*sd).ghost_size;
    int id_recv_up = ((*sd).ghost_size + ((*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2]));
    int id_send_down = (((*sd).grid_l[1] * (*sd).grid_l[2]) * ((*sd).grid_l[0] - 1));
    int id_recv_down = 0;

    MPI_Status status;

    if (type == FDM)
    {
        // send to top slice of a subdomain to the subdomain above it
        MPI_Sendrecv(&((*sd).u_next[id_send_up]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT, (*sd).neighbors[Down], Up, 
                    &((*sd).u_next[id_recv_up]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT, (*sd).neighbors[Up], Up, 
                    sd->COMM_FDM, &status);

        // send the bottom slice of a subdomain to the subdomain below it
        MPI_Sendrecv(&((*sd).u_now[id_send_down]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT, (*sd).neighbors[Up], Down, 
                    &((*sd).u_now[id_recv_down]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_FLOAT, (*sd).neighbors[Down], Down, 
                    sd->COMM_FDM, &status);
    }
    else if (type == Shape)
    {
        // send to top ghost of a subdomain to the subdomain above it
        MPI_Sendrecv(&((*sd).shape_next[id_send_up]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, (*sd).neighbors[Down], Up, 
                    &((*sd).shape_next[id_recv_up]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, (*sd).neighbors[Up], Up, 
                    sd->COMM_FDM, &status);

        // send the bottom ghost of a subdomain to the subdomain below it
        MPI_Sendrecv(&((*sd).shape_now[id_send_down]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, (*sd).neighbors[Up], Down, 
                    &((*sd).shape_now[id_recv_down]), (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, (*sd).neighbors[Down], Down, 
                    sd->COMM_FDM, &status);
    }
} // end void ShareGhosts(Subdomain *sd, int type)

void DetermineSubdomainGridBounds(Subdomain *sd)
{
    // assign the process for this subdomain cartesian coordinates.
    if (sd->n_dims == 2)
    {
        sd->coords[0] = floor(sd->rank / sd->n_proc_dim[1]);
        sd->coords[1] = sd->rank - sd->coords[0] * sd->n_proc_dim[1];
        sd->coords[2] = 0;
    }
    else if (sd->n_dims == 3)
    {
        // to work out these equations, think of converting 10,000 seconds to hr:min:s (z:y:x)
        sd->coords[2] = floor(sd->rank / (sd->n_proc_dim[0] * sd->n_proc_dim[1]));
        sd->coords[1] = floor((sd->rank - sd->coords[2] * sd->n_proc_dim[0] * sd->n_proc_dim[1]) / sd->n_proc_dim[0]);
        sd->coords[0] = (sd->rank - sd->coords[2] * sd->n_proc_dim[0] * sd->n_proc_dim[1]) - sd->n_proc_dim[0] * sd->coords[1];
    }

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
            sd->bounds_l[2 * n + 3] = 1;
        }
    }
} // end void DetermineSubdomainGridBounds(Subdomain *sd)

void CollectSubdomainData(Subdomain *sd)
{       
    int offset = sd->ghost_size;
    int global_id = ((sd->bounds_l[0] * sd->grid_g[1] + sd->bounds_l[2]) * sd->grid_g[2] + sd->bounds_l[4]);  
    //MPI_Gather(&((*sd).shape_next[offset]), (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, 
      //         &((*sd).shape_g[global_id]), (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT,
        //       ROOT, MPI_COMM_WORLD);

    MPI_Gather(&((*sd).shape_now[offset]), (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT, 
               &((*sd).shape_g[global_id]), (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT,
               ROOT, MPI_COMM_WORLD);
} // end void CollectSubdomainData(Subdomain *subdomain)

void AllocateArraysFDM(Subdomain *sd)
{
    // size of memory to allocate
    int local_size, global_size;

    if (sd->n_proc == 1) { sd->ghost_size = 0; } // sequential, so no ghost cells
    else if (sd->n_proc !=1 && sd->n_dims == 2) { sd->ghost_size = sd->grid_l[1]; } // ghost cell is size of row/column
    else if (sd->n_proc !=1 && sd->n_dims == 3) { sd->ghost_size = sd->grid_l[1] * sd->grid_l[2]; } // ghost cell is size of plane

    local_size = ((*sd).grid_l[0]) * ((*sd).grid_l[1]) * ((*sd).grid_l[2]) + (2 * sd->ghost_size);
    global_size = ((*sd).grid_g[0]) * ((*sd).grid_g[1]) * ((*sd).grid_g[2]);

    // local shape arrays
    (*sd).shape_now = (int*)malloc(sizeof(int) * local_size);
    (*sd).shape_next = (int*)malloc(sizeof(int) * local_size);
    memset((*sd).shape_now, 0, sizeof(int) * local_size);
    memset((*sd).shape_next, 0, sizeof(int) * local_size);

    // local solution arrays
    (*sd).u_next = (float*)malloc(sizeof(float) * local_size);
    (*sd).u_now = (float*)malloc(sizeof(float) * local_size);
    memset((*sd).u_next, 0, sizeof(float) * local_size);
    memset((*sd).u_now, 0, sizeof(float) * local_size);

    // global shape array
    (*sd).shape_g = (int*)malloc(sizeof(int) * global_size);
    memset((*sd).shape_g, 0, sizeof(int) * global_size);
    
    (*sd).u_global = (float*)malloc(sizeof(float) * global_size);
    memset((*sd).u_global, 0, sizeof(float) * global_size);
} // end void AllocateArraysFDM(Subdomain sd)

void SubdomainCleanUp(Subdomain *sd)
{
    free((*sd).shape_now);
    free((*sd).shape_next);
    free((*sd).u_now);
    free((*sd).u_next);

    if (sd->rank == ROOT)
    {
        free((*sd).shape_g);
        free((*sd).u_global);
    }

    free(sd);
} // void SubdomainCleanUp(Subdomain *sd)

void CoordShift(Subdomain *sd, float r)
{
    //fprintf(stdout, "coord shift\n");
    int local_id; // the computations below are shifted to global coordinates, this converts back to local index
    int offset = sd->ghost_size;
    fprintf(stdout, "offset is %d\n", offset);
    float r_shifted; // squared length of shifted position vector of a fdm grid_g cell
    for (int i = sd->bounds_l[0]; i < sd->bounds_l[1]; i ++)
    {
        for (int j = sd->bounds_l[2]; j < sd->bounds_l[3]; j++)
        {
            for (int k = sd->bounds_l[4]; k < sd->bounds_l[5]; k++)
            {
                local_id = ID(sd, (i - sd->bounds_l[0]), (j - sd->bounds_l[2]), (k - sd->bounds_l[4])) + offset;
                
                r_shifted = pow(i - (sd->grid_g[0] / 2), 2.0) + pow(j - (sd->grid_g[1] / 2), 2.0);
                if (sd->n_dims == 3) { r_shifted += pow(k - (sd->grid_g[2] / 2), 2.0); }
                
                if (r_shifted <= r*r) // inside of circle
                {
                    fprintf(stdout, "rank %d: r_shifted = %f\n", sd->rank, r_shifted);
                    fprintf(stdout, "rank %d: sd->n_proc_dim[0] = %d\n", sd->rank, sd->n_proc_dim[0]);
                    fprintf(stdout, "rank %d: sd->n_proc_dim[1] = %d\n", sd->rank, sd->n_proc_dim[1]);
                    fprintf(stdout, "coordinates (%d, %d, %d)\n", i, j, k);
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
            } // end j loop
        } // end i loop
    }
    MPI_Barrier(MPI_COMM_WORLD);
} // end void CoordShift(Subdomain sd, float r)