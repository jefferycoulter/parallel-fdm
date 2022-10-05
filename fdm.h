#ifndef FDM_INCL
#define FDM_INCL

/**
 * @brief store space/time StepSize in a single struct to reduce parameter passing
 */
typedef struct 
{
    float dx, dy, dz; // spatial step size
    float dt; // temporal step size
} StepSize;

/**
 * @brief store data arrays in a single struct to reduce parameter passing
 */
typedef struct
{
    float *next; // stores data after FDM iteration
    float *now; // used to compute next iteration
} Data;

/**
 * @brief specifies number of dimensions and the size of the global simulation 
 * domain along each direction
 * @param n
 */
typedef struct
{
    int n; // number of dimensions
    int x, y, z; // dimension sizes
} Dimension;

/**
 * @brief run finite difference method for a given number of iterations
 * @param next result of FDM computation after the next iteration
 * @param data_current either initial data, or data from currnt iteration
 * @param dt time step size
 * @param time_steps number of time steps
 * @param dx space step size
 * @param cells number of cells along each axis
 */
void ComputeFD(Data data, StepSize dX, int time_steps, int cells[3]);

/**
 * @brief compute a single temporal iteration of forward-time central-space (FTCS) finite difference 
 * scheme for a parabolic PDE (2D diffusion equation right now)
 * @param next result of FDM computation after the next iteration
 * @param data_current either initial data, or data from currnt iteration
 * @param dt time step size
 * @param dx space step size
 * @param dim number of spatial dimensions
 */
void FTCS(Data data, StepSize dX, int cells[3]);

/**
 * @brief apply boundary conditions to the domain
 * @param data simulation domain where boundary conditions will be applied
 * @param cells spatial discretizaion along each axis
 */
void SetBoundaryConditions(Data data, int cells[3]);

/**
 * @brief apply initial conditions to the domain
 * @param data simulation domain where boundary conditions will be applied
 * @param cells spatial discretizaion along each axis
 */
void SetInitialConditions(float *data, int cells[3]);

/**
 * @brief discretize the process's domain (a subdomain)
 * @param dims size of the domain in the 
 * @param dx spatial step size
 * @param cells resulting array returns number of discretizations along each direction
 */
void DiscretizeSubdomain(Dimension dim, StepSize dX, int cells[3]);

/**
 * @brief allocate memory for domain arrays
 * @param next the "next" iteration to solve
 * @param data_current the "current" iteration solution for computing next iteration
 * @param cells spatial discretization along each axis
 */
void GenerateSpatialArrays(Data *data, int cells[3]);

#endif // FDM_INCL