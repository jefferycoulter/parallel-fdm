#ifndef PARALLEL_INCL
#define PARALLEL_INCL

#include <mpi.h>

#define ROOT 0 // root process

/**
 * @brief struct containing items related to subdomain construction and usage
 */
typedef struct
{
    int n_proc; // total number of processors
    int n_proc_dim[3]; // number of processors along each dimension
    int rank; // rank of process that owns the subdomain in MPI_COMM_WORLD
    MPI_Datatype subdomain_type;
    int coords[3];
    int *send_counts;
    int *displs;

    int n_dims; // number of dimensions
    int dims_g[3]; // global domain length along each dimension
    int dims_l[3]; // local/subdomain length along each dimension

    int bounds_l[6]; // bounds of subdomain (local bounds) filled like [i_start, i_end, j_start, j_end, k_start, k_end]
    int grid_g[3]; // number of global fdm spatial grid cells along each dimension
    int grid_l[3]; // number of local fdm spatial grid cells along each dimension

    float dx, dy, dz; // spatial step size
    float dt; // temporal step size
    float mu_x, mu_y, mu_z; // stability factor, mu_i = dt / di**2 where i = x, y, z

    int *shape_now; // shape array local to a process. used for determining fdm boundaries in case of irregular geometry
    int *shape_next;
    int *shape_g; // global shape array 
    float *u_now; // current result, used to compute next fdm iteration
    float *u_next; // stores data after fdm iteration
    float *u_global; // global solution array
} Subdomain;

/**
 * @brief create a subdomain structure on each process.  this handles member initializations,
 * memory allocations, and distributes processors amongst the spatial dimensions
 * @param nproc number of processors being used
 * @param rank rank of the processor that the subdomain belongs to
 * @return Subdomain* 
 */
Subdomain *CreateSubdomain(int nproc, int rank);

/**
 * @brief get input parameters to create the subdomain. this is probably not going to last. can initialize
 * the subdomain in the main function with the Input struct defined in io.h
 * @param subdomain subdomain to supply input to
 * @param nproc number of processors being used
 * @param rank rank of the processor that the subdomain belongs to
 */
void GetInput(Subdomain *subdomain, int nproc, int rank);

/**
 * @brief divide processors amongst the spatial dimensions in an optimal way
 * @param subdomain 
 */
void SplitProcessorsAlongDims(Subdomain *subdomain);

/**
 * @brief determine the subdomain upper and lower bounds relative to the global grid
 * @param subdomain 
 */
void DetermineSubdomainGridBounds(Subdomain *subdomain);

/**
 * @brief allocate memory for FDM arrays
 * @param subdomain the subdomain to generate arrays for
 */
void AllocateArraysFDM(Subdomain *subdomain);

/**
 * @brief create a custom MPI_Datatype for subdomains. this is used to ensure that when data from each
 * subdomain is sent to ROOT for writing results to a file, the subdomains are placed in the correct order
 * @param subdomain 
 */
void CreateSubdomainType(Subdomain *subdomain);

void CreateSubdomainType2(Subdomain *sd);

void SetupCollectSubdomainData(Subdomain *subdomain);

/**
 * @brief essentially just a wrapper for MPI_Gather right now. will add more later
 * @param subdomain subdomain to collect from
 */
void CollectSubdomainData(Subdomain *subdomain);

void SubdomainCleanUp(Subdomain *sd);

/**
 * @brief shift the coordinates in a given subdomain so that the center of the global subdomain (dim_X / 2, dim_y / 2) 
 * is moved to the origin (0,0).  this is used for determining boundaries in the case of irregular geometry
 * @param subdomain subdomain to shift
 * @param radius the radius of the irregular geometry
 */
void CoordShift(Subdomain *subdomain, float radius);

#endif // PARALLEL_INCL