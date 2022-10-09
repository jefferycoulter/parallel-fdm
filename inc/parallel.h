#ifndef PARALLEL_INCL
#define PARALLEL_INCL

#define ROOT 0 // root process

/**
 * @brief struct containing items related to subdomain construction and usage
 */
typedef struct
{
    int n_proc; // total number of processors
    int n_proc_dim[3]; // number of processors along each dimension
    int rank; // rank of process that owns the subdomain

    int n_dims; // number of dimensions
    int dims_g[3]; // global domain length along each dimension
    int dims_l[3]; // local/subdomain length along each dimension

    int bounds_l[6]; // bounds of subdomain (local bounds), i.e. i_start, i_end, j_start, etc.
    int *shape_arr_l; // shape array local to a process. used for determining fdm boundaries in case of irregular geometry
    int *shape_arr_g; // global shape array 
    int grid_g[3]; // number of global fdm spatial grid cells along each dimension
    int grid_l[3]; // number of local fdm spatial grid cells along each dimension

    float dx, dy, dz; // spatial step size
    float dt; // temporal step size
    float mu_x, mu_y, mu_z; // stability factor, mu_i = dt / di**2 where i = x, y, z

    float *u_now; // current result, used to compute next fdm iteration
    float *u_next; // stores data after fdm iteration
    float *u_global; // global solution array
} Subdomain;

void PrepareSubdomains(Subdomain *subdomain);

void SplitProcessorsAlongDims(Subdomain *subdomain);

void SetTopology(Subdomain subdomain);

/**
 * @brief get the boundaries for each subdomain, i.e. x_start, x_end, y_start, etc.
 * @param subdomain the subdomain belonging to a given process
 */
void GetSubdomainGridBounds(Subdomain *subdomain);

/**
 * @brief compute the local size of the subdomain along each dimension
 * @param subdomain 
 */
void ComputeSubdomainGrid(Subdomain *subdomain);

/**
 * @brief shift the coordinates in a given subdomain so that the center of the global subdomain (dim_X / 2, dim_y / 2) 
 * is moved to the origin (0,0).  this is used for determining boundaries in the case of irregular geometry
 * @param subdomain subdomain to shift
 * @param radius the radius of the irregular geometry
 */
void CoordShift(Subdomain *subdomain, float radius);

#endif // PARALLEL_INCL