#ifndef FDM_INCL
#define FDM_INCL

#include "parallel.h"

/**
 * @brief convert point in 3D space to a corresponding index in linear memory
 * @param sd subdomain containing the point
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define ID(sd, i, j, k) (i * (*sd).grid_l[1] + j) * (*sd).grid_l[2] + k

/**
 * @brief boundary conditions that can be used for fdm. outside region of interest (i.e. no FDM computations here),
 * inside region of interest (i.e. where FDM computations occur), or boundary of region of interestion
 */
enum BC{Dirichlet, VonNeumann};

/**
 * @brief location of fdm grid cell
 */
enum Location{Outside, Inside, Boundary};

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
 * @param id index of grid point
 */
void InteriorFD(Subdomain *subdomain, int id);

/**
 * @brief compute boundary terms depending on whether dirichlet or von neuman is being used
 * @param subdomain the subdomain on which to compute finite differences
 * @param id index of grid point
 */
void BoundaryFD(Subdomain *subdomain, int bc, int id);

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
 * @brief determine the step sizes along each dimension from global size and the
 * number of finite difference grid_g cells
 * @param sd 
 */
void DetermineStepSizes(Subdomain *sd);

/**
 * @brief finite difference in x direction
 * @param sd subdomain
 * @param id index of grid point
 */
#define FDX(sd, u, id)     (1 - 2 * (*sd).mu_x) * u[id]                \
                        + (*sd).mu_x * (u[id + (*sd).grid_l[1]])    \
                        + (*sd).mu_x * (u[id - (*sd).grid_l[1]])

/**
 * @brief finite difference in y direction
 * @param sd subdomain
 * @param id index of grid point
 */
#define FDY(sd, u, id)     - 2 * (*sd).mu_y * u[id]      \
                        + (*sd).mu_y * (u[id + 1])    \
                        + (*sd).mu_y * (u[id - 1])

/**
 * @brief finite difference in z direction
 * @param sd subdomain
 * @param id index of grid point
 */
#define FDZ(sd, u, id)     - 2 * (*sd).mu_z * u[id]                                      \
                        + (*sd).mu_z * (u[id + ((*sd).grid_l[1] * (*sd).grid_l[2])])  \
                        + (*sd).mu_z * (u[id - ((*sd).grid_l[1] * (*sd).grid_l[2])])

/**
 * @brief apply finite difference at point (i, j, k)
 * @param sd subdomain
 * @param id index of grid point
 */
#define FD(sd, u, id) ((sd->n_dims) == 2 ? FDX(sd, u, id) + FDY(sd, u, id) : FDX(sd, u, id) + FDY(sd, u, id) + FDZ(sd, u, id))

#endif // FDM_INCL