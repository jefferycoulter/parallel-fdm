#ifndef FDM_INCL
#define FDM_INCL

#include "parallel.h"
#include <math.h>

/**
 * @brief boundary conditions that can be used for fdm. outside region of interest (i.e. no FDM computations here),
 * inside region of interest (i.e. where FDM computations occur), or boundary of region of interestion
 */
enum BC{Dirichlet, VonNeumann};

/**
 * @brief convert point in 3D space to a corresponding index in linear memory
 * @param sd subdomain containing the point
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define ID(sd, i, j, k) ((i * sd->grid_l[1]) + j) * sd->grid_l[2] + k

/**
 * @brief location of fdm grid cell
 * 
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

#define FDX(sd, i, j, k, off)   (1 - 2 * sd->mu_x) * sd->u_now[off + (i * sd->grid_l[1] + j) * sd->grid_l[2] + k]                             \
                                + sd->mu_x * (sd->u_now[off + (i * sd->grid_l[1] + j + sd->grid_l[0]) * sd->grid_l[2] + k])                   \
                                + sd->mu_x * (sd->u_now[off + (i * sd->grid_l[1] + j - sd->grid_l[0]) * sd->grid_l[2] + k])

#define FDY(sd, i, j, k, off)   (1 - 2 * sd->mu_y) * sd->u_now[off + (i * sd->grid_l[1] + j) * sd->grid_l[2] + k]                             \
                                + sd->mu_y * (sd->u_now[off + (i * sd->grid_l[1] + j + 1) * sd->grid_l[2] + k])                               \
                                + sd->mu_y * (sd->u_now[off + (i * sd->grid_l[1] + j - 1) * sd->grid_l[2] + k])

#define FDZ(sd, i, j, k, off)   (1 - 2 * sd->mu_z) * sd->u_now[off + (i * sd->grid_l[1] + j) * sd->grid_l[2] + k]                             \
                                + sd->mu_z * (off + sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k + sd->grid_l[0] * sd->grid_l[1]])   \
                                + sd->mu_z * (off + sd->u_now[(i * sd->grid_l[1] + j) * sd->grid_l[2] + k - sd->grid_l[0] * sd->grid_l[1]])

#define FD(sd, i, j, k, off) ((k) == 0 ? FDX(sd, i, j, k, off) + FDY(sd, i, j, k, off) : FDX(sd, i, j, k, off) + FDY(sd, i, j, k, off) + FDZ(sd, i, j, k, off))

#endif // FDM_INCL