#ifndef FDM_INCL
#define FDM_INCL

#include "parallel.h"

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
 * @param r row index
 * @param c column index
 * @param d depth index
 */
void InteriorFD(Subdomain *subdomain, int r, int c, int d);

/**
 * @brief compute boundary terms depending on whether dirichlet or von neuman is being used
 * @param subdomain the subdomain on which to compute finite differences
 * @param r row index
 * @param c column index
 * @param d depth index
 */
void BoundaryFD(Subdomain *subdomain, int bc, int r, int c, int d);

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

/**
 * @brief macro to decide if the point is interior of the domain or a boundary point
 * @param sd subdomain
 * @param bc boundary conditions
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define InteriorOrBoundary(sd, bc, r, c, d) switch (sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c])    \
                                            {                                                                       \
                                            case OUTSIDE:                                                           \
                                                break;                                                              \
                                            case INSIDE:                                                            \
                                                InteriorFD(sd, c, r, d);                                            \
                                                break;                                                              \
                                            case BOUNDARY:                                                          \
                                                BoundaryFD(sd, bc, c, r, d);                                        \
                                                break;                                                              \
                                            } 

/**
 * @brief finite difference in x direction
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define FDX(sd, r, c, d) (1 - 2 * sd->mu_x) * sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c]    \
                         + sd->mu_x * (sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c + 1]       \
                                      + sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c - 1])

/**
 * @brief finite difference in y direction
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define FDY(sd, r, c, d) (1 - 2 * sd->mu_y) * sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c]            \
                         + sd->mu_y * (sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c + sd->grid_l[0]]   \
                                      + sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c - sd->grid_l[0]])

/**
 * @brief finite difference in z direction
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define FDZ(sd, r, c, d) (1 - 2 * sd->mu_z) * sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c]                                \
                         + sd->mu_z * (sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c + sd->grid_l[0] * sd->grid_l[1]]       \
                                      + sd->u_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c - sd->grid_l[0] * sd->grid_l[1]])

/**
 * @brief increment the rows, columns, and depth if necessary
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define Increment(i, r, c, d)  c += 1; /* increment columns */                                              \
                            if ((i % sd->grid_l[0]) == 0) /* reached end of a row */                        \
                            {                                                                               \
                                r += 1; /* increment row */                                                 \
                                c = 0;  /* reset columns */                                                 \
                            }                                                                               \
                            if (sd->n_dims == 3)                                                            \
                            {                                                                               \
                                if ((i % sd->grid_l[0] * sd->grid_l[1]) == 0) /* reached end of plane */    \
                                {                                                                           \
                                    d += 1; /* increment depth */                                           \
                                    c = 0;  /* reset columns */                                             \
                                    r = 0;  /* reset rows */                                                \
                                }                                                                           \
                            }

#endif // FDM_INCL