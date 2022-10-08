#ifndef PARALLEL_INCL
#define PARALLEL_INCL

/**
 * @brief struct containing items related to subdomain construction and usage
 */
typedef struct
{
    int n_proc; // total number of processors
    int rank; // rank of process that owns the subdomain

    int n_dims; // number of dimensions
    int dims_g[3]; // global domain length along each dimension
    int dims_l[3]; // local/subdomain length along each dimension

    int bounds_l[6]; // bounds of subdomain (local bounds), i.e. x_start, x_end, y_start, etc.
    int *shape_arr; // shape array used for determining fdm boundaries in case of irregular geometry
    int grid[3]; // fdm spatial grid

    float dx, dy, dz; // spatial step size
    float dt; // temporal step size
    float mu_x, mu_y, mu_z; // stability factor, mu_i = dt / di**2 where i = x, y, z

    float *u_now; // current result, used to compute next fdm iteration
    float *u_next; // stores data after fdm iteration
} Subdomain;

/**
 * @brief get the boundaries for each subdomain, i.e. x_start, x_end, y_start, etc.
 * @param subdomain the subdomain belonging to a given process
 */
void GetSubdomainBounds(Subdomain subdomain);


/**
 * @brief compute the local size of the subdomain along each dimension
 * @param subdomain 
 */
void ComputeSubdomainDims(Subdomain subdomain);

/**
 * @brief shift the coordinates in a given subdomain so that the center of the global subdomain is at the
 * global origin.  this is used for determining boundaries in the case of irregular geometry
 * @param subdomain 
 */
void CoordShift(Subdomain subdomain, float r);

#endif // PARALLEL_INCL