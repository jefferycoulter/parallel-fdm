#ifndef SHAPE_INCL
#define SHAPE_INCL

#include "fdm.h"

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
 * @brief shift the coordinates in a given subdomain so that the center of the global subdomain (dim_X / 2, dim_y / 2) 
 * is moved to the origin (0,0).  this is used for determining boundaries in the case of irregular geometry
 * @param subdomain subdomain to shift
 * @param radius the radius of the irregular geometry
 */
void CoordShift(Subdomain *subdomain, float radius);

/**
 * @brief laplace filter in x direction
 * @param sd subdomain
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define LX(sd, id)  -2 * (*sd).shape_now[id]                   \
                    + (*sd).shape_now[id + (*sd).grid_l[1]]    \
                    + (*sd).shape_now[id - (*sd).grid_l[1]]

/**
 * @brief laplace filter in y direction
 * @param sd subdomain
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define LY(sd, id)  -2  * (*sd).shape_now[id]   \
                    + (*sd).shape_now[id + 1]   \
                    + (*sd).shape_now[id - 1]

/**
 * @brief laplace filter in z direction
 * @param sd subdomain
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define LZ(sd, id)  -2 * (*sd).shape_now[id]                                        \
                    + (*sd).shape_now[id + ((*sd).grid_l[1] * (*sd).grid_l[2])]     \
                    + (*sd).shape_now[id - ((*sd).grid_l[1] * (*sd).grid_l[2])] 

/**
 * @brief apply the laplace filter at point (i, j, k)
 * @param sd subdomain
 * @param i x index
 * @param j y index
 * @param k z index
 * @param off offset due to ghost cell
 */
#define Laplace(sd, id, k) ((k) == 0 ? LX(sd, id) + LY(sd, id) : LX(sd, id) + LY(sd, id) + LZ(sd, id))

/**
 * @brief after performing the laplace operation, assign values to the boundary, the interior domain, and the exterior domain.
 * @param sd subdomain
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define AssignValue(sd, id) switch (sd->shape_next[id])                                                                       \
                            {                                                                                                 \
                                case 0: /* if the new value is zero, then the previous value was either INSIDE or OUTSIDE*/   \
                                    switch (sd->shape_now[id])                                                                \
                                    {                                                                                         \
                                        case Outside: /* if previous value was OUTSIDE, then new value is still outside */    \
                                            sd->shape_next[id] = Outside;                                                     \
                                            break;                                                                            \
                                                                                                                              \
                                        case Inside: /* if previous value was INSIDE, then switch it back to inside */        \
                                            sd->shape_next[id] = Inside;                                                      \
                                            break;                                                                            \
                                    }                                                                                         \
                                    break;                                                                                    \
                                default: /* if the new value is not zero, then it is a boundary */                            \
                                    sd->shape_next[id] = Boundary;                                                            \
                                    break;                                                                                    \
                            }

#endif // SHAPE_INCL