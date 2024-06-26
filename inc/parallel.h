#ifndef PARALLEL_INCL
#define PARALLEL_INCL

#include <mpi.h>

#define ROOT 0 // root process

/**
 * @brief enumeration for finding neighbor processors
 */
enum Neighbors{Up, Down};

/**
 * @brief struct containing items related to subdomain construction and usage
 */
typedef struct
{
    int n_proc; // total number of processors
    int n_proc_dim[3]; // number of processors along each dimension
    int rank; // rank of process that owns the subdomain in MPI_COMM_WORLD
    MPI_Comm COMM_FDM; // finite difference cartesian communicator
    int coords[3]; // cartesian coordinates of process
    int neighbors[2]; // neighbor processes

    int n_dims; // number of dimensions
    int dims_g[3]; // global domain length along each dimension

    int bounds_l[6]; // bounds of subdomain (local bounds), i.e. i_start, i_end, j_start, etc.
    int grid_g[3]; // number of global fdm spatial grid cells along each dimension
    int grid_l[3]; // number of local fdm spatial grid cells along each dimension

    float dx, dy, dz; // spatial step size
    float dt; // temporal step size

    float mu_x, mu_y, mu_z; // stability factor, mu_i = dt / di**2 where i = x, y, z

    int ghost_size; // size in memory of ghost cells

    int *shape_now; // shape array local to a process. used for determining fdm boundaries in case of irregular geometry
    int *shape_next;
    int *shape_g; // global shape array -- this might not be necessary. just used for debugging i think
    
    // species 1
    float k1;
    float *u_now; // current result, used to compute next fdm iteration
    float *u_next; // stores data after fdm iteration
    float *u_global; // global solution array

    // species 2
    float k2;
    float *v_now; // current result, used to compute next fdm iteration
    float *v_next; // stores data after fdm iteration
    float *v_global; // global solution array

    // species 3
    float k3;
    float *uv_now; // current result, used to compute next fdm iteration
    float *uv_next; // stores data after fdm iteration
    float *uv_global; // global solution array
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
 * @brief get input parameters to create the subdomain
 * @param subdomain subdomain to supply input to
 * @param nproc number of processors being used
 * @param rank rank of the processor that the subdomain belongs to
 */
void LoadParameters(Subdomain *subdomain, int nproc, int rank);

/**
 * @brief if the given global dimensions of the domain aren't divisible by the number of processors
 * then slightly adjust them
 * @param subdomain
 * @param dim dimension which should be resized 
 */
void ResizeDomain(Subdomain *subdomain, int dim);

/**
 * @brief divide processors amongst the spatial dimensions in an optimal way
 * @param subdomain 
 */
void SplitProcessorsAlongDims(Subdomain *subdomain);

/**
 * @brief get the neighbor processors for communication
 * @param subdomain the subdomain to get neighbors for
 */
void GetNeighbors(Subdomain *subdomain);

/**
 * @brief send and recieve ghost cells from neighboring processors
 * @param subdomain subdomain to send/recieve
 * @param type FDM or Shape. specifies which arrays should share ghost cells
 */
void ShareGhosts(Subdomain *subdomain, int type);

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
 * @brief collect local data to root process for writing to a file
 * @param subdomain the subdomain to collect from
 * @param type FDM or Shape. specifies which arrays should be collected
 * @param time time step at which data is collected
 */
void CollectSubdomainData(Subdomain *subdomain, int type, int time);

/**
 * @brief free memory on the process
 * @param subdomain the subdomain to clean up
 */
void SubdomainCleanUp(Subdomain *subdomain);

#endif // PARALLEL_INCL