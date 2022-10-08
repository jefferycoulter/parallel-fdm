#ifndef FDM_INCL
#define FDM_INCL

#include "parallel.h"

#include <math.h>

/**
 * @brief boundary conditions that can be used for fdm
 */
enum BC{Dirichlet, VonNeumann};

/**
 * @brief run finite difference method for a given number of iterations
 * @param next result of FDM computation after the next iteration
 * @param data_current either initial data, or data from currnt iteration
 * @param dt time step size
 * @param time_steps number of time steps
 * @param dx space step size
 * @param cells number of cells along each axis
 */
void ComputeFD(Subdomain subdomain, int bc, int time_steps);

/**
 * @brief compute a single temporal iteration of forward-time central-space (FTCS) finite difference 
 * scheme for a parabolic PDE (2D diffusion equation right now)
 * @param next result of FDM computation after the next iteration
 * @param data_current either initial data, or data from currnt iteration
 * @param dt time step size
 * @param dx space step size
 * @param dim number of spatial dimensions
 */
void FTCS(Subdomain subdomain, int bc);

/**
 * @brief compute interior points. same for both dirichlet and von neumann
 * @param data data to compute finite differences on
 * @param dX spatial and temporal step sizes
 * @param cells number of cells along each dimension
 */
void InteriorFD(Subdomain subdomain);

/**
 * @brief compute boundary terms depending on whether dirichlet or von neuman is being used
 * @param data data to compute finite differences on
 * @param bc 
 * @param dX spatial and temporal step sizes
 * @param cells number of cells along each dimension
 */
void BoundaryFD(Subdomain subdomain, int bc);

void CreateShapeArray(Subdomain subdomain, float radius);

/**
 * @brief apply boundary conditions to the domain
 * @param data simulation domain where boundary conditions will be applied
 * @param cells spatial discretizaion along each axis
 */
void SetBoundaryConditions(Subdomain subdomain);

/**
 * @brief apply initial conditions to the domain
 * @param data simulation domain where boundary conditions will be applied
 * @param cells spatial discretizaion along each axis
 */
void SetInitialConditions(Subdomain subdomain);

/**
 * @brief discretize the process's domain (a subdomain)
 * @param dims size of the domain in the 
 * @param dx spatial step size
 * @param cells resulting array returns number of discretizations along each direction
 */
void DiscretizeSubdomain(Subdomain subdomain);

/**
 * @brief allocate memory for domain arrays
 * @param next the "next" iteration to solve
 * @param data_current the "current" iteration solution for computing next iteration
 * @param cells spatial discretization along each axis
 */
void GenerateArraysFDM(Subdomain subdomain);

#endif // FDM_INCL