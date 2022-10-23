#ifndef PARALLEL_INCL
#define PARALLEL_INCL

#include <mpi.h>

#define ROOT 0 // root process

/**
 * @brief enumeration for finding neighbor processors
 */
enum Neighbors{Front, Back, Left, Right, Up, Down};

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

    int bounds_l[6]; // bounds of subdomain (local bounds), i.e. i_start, i_end, j_start, etc.
    int grid_g[3]; // number of global fdm spatial grid cells along each dimension
    int grid_l[3]; // number of local fdm spatial grid cells along each dimension

    float dx, dy, dz; // spatial step size
    float dt; // temporal step size
    float mu_x, mu_y, mu_z; // stability factor, mu_i = dt / di**2 where i = x, y, z
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

void CreateSubdomainType(Subdomain *subdomain);

void CreateSubdomainType2(Subdomain *subdomain);

void SetupCollectSubdomainData(Subdomain *subdomain);

/**
 * @brief gather data from all processes to root
 */
#define CollectSubdomainData(sd, loc, glo)                                                      \
    MPI_Gather(loc, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT,  \
                glo, (*sd).grid_l[0] * (*sd).grid_l[1] * (*sd).grid_l[2], MPI_INT,    \
                ROOT, MPI_COMM_WORLD);                                                          \

//void CollectSubdomainData(Subdomain *subdomain);

void SubdomainCleanUp(Subdomain *subdomain);


#endif // PARALLEL_INCL