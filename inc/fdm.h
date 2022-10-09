#ifndef FDM_INCL
#define FDM_INCL

#include "parallel.h"
#include <math.h>

#define OUTSIDE 0 // outside region of interest (i.e. no FDM computations here)
#define INSIDE 1 // region of interest (i.e. where FDM computations occur)
#define BOUNDARY 2 // boundary of region of interestion

/**
 * @brief boundary conditions that can be used for fdm
 */
enum BC{Dirichlet, VonNeumann};

/**
 * @brief run finite difference method for a given number of iterations
 * @param subdomain the subdomain on which to compute finite differences
 * @param bc boundary conditions (Dirichlet or VonNeumann)
 * @param time_steps max number of iterations for computations
 */
void ComputeFD(Subdomain *subdomain, int bc, int time_steps);

/**
 * @brief compute a single temporal iteration of forward-time central-space (FTCS) finite difference 
 * scheme for a parabolic PDE (2D diffusion equation right now)
 * @param subdomain the subdomain on which to compute finite differences
 * @param bc boundary conditions (Dirichlet or VonNeumann)
 */
void FTCS(Subdomain *subdomain, int bc);

/**
 * @brief compute interior points. same for both dirichlet and von neumann
 * @param subdomain the subdomain on which to compute finite differences
 * @param i x index
 * @param j y index
 * @param z index (zero if 2D)
 */
void InteriorFD(Subdomain *subdomain, int i, int j, int k);

/**
 * @brief compute boundary terms depending on whether dirichlet or von neuman is being used
 * @param subdomain the subdomain on which to compute finite differences
 * @param i x index
 * @param j y index
 * @param z index (zero if 2D)
 */
void BoundaryFD(Subdomain *subdomain, int bc, int i, int j, int k);

/**
 * @brief create a shape array.  in the case of irregular geometry, this corresponds to which locations
 * are inside the area of interest (i.e. where FDM computation occurs), which locations are boundaries,
 * and which locations are outside the area of interest
 * @param subdomain subdomain to create shape array in
 * @param radius radius of domain of interest
 */
void CreateShapeArray(Subdomain *subdomain, float radius);

/**
 * @brief apply a laplace filter (discrete laplace operator) to detect the boundary of the area
 * of interest (i.e. where FDM computation occurs)
 * @param subdomain subdomain to apply filter on
 */
void ApplyLaplaceFilter(Subdomain *subdomain);

/**
 * @brief apply boundary conditions to the domain
 * @param subdomain the subdomain in which to set boundary conditions
 */
void SetBoundaryConditions(Subdomain *subdomain);

/**
 * @brief apply initial conditions to the domain
 * @param subdomain the subdomain in which to set initial conditions
 */
void SetInitialConditions(Subdomain *subdomain);

/**
 * @brief discretize the process's domain (a subdomain)
 * @param subdomain the subdomain to discretize
 */
void DiscretizeSubdomain(Subdomain *subdomain);

/**
 * @brief determine the step sizes along each dimension from global size and the
 * number of finite difference grid_g cells
 * @param sd 
 */
void DetermineStepSizes(Subdomain *sd);

#endif // FDM_INCL