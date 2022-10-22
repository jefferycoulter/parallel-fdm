#ifndef SHAPE_INCL
#define SHAPE_INCL

#include "parallel.h"

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
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define LX(sd, r, c, d) -2  * sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c]  \
                        + sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c + 1]  \
                        + sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c - 1]

/**
 * @brief laplace filter in y direction
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define LY(sd, r, c, d) -2 * sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c]               \
                        + sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c + sd->grid_l[0]]  \
                        + sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c - sd->grid_l[0]]

/**
 * @brief laplace filter in z direction
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define LZ(sd, r, c, d) -2 * sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c]                               \
                        + sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c + sd->grid_l[0] * sd->grid_l[1]]  \
                        + sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c - sd->grid_l[0] * sd->grid_l[1]]

/**
 * @brief apply the laplace filter at point (r, c, d)
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define Laplace(sd, r, c, d) ((d) == 0 ? LX(sd, r, c, d) + LY(sd, r, c, d) : LX(sd, r, c, d) + LY(sd, r, c, d) + LZ(sd, r, c, d))

/**
 * @brief after performing the laplace operation, assign values to the boundary, the interior domain, and the exterior domain.
 * @param sd subdomain
 * @param r row index
 * @param c column index
 * @param d depth index
 */
#define AssignValue(sd, r, c, d)    switch (sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c])                              \
                                    {                                                                                                 \
                                        case 0: /* if the new value is zero, then the previous value was either INSIDE or OUTSIDE*/   \
                                            switch (sd->shape_now[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c])                       \
                                            {                                                                                         \
                                                case OUTSIDE: /* if previous value was OUTSIDE, then new value is still outside */    \
                                                    sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = OUTSIDE;            \
                                                    break;                                                                            \
                                                                                                                                      \
                                                case INSIDE: /* if previous value was INSIDE, then switch it back to inside */        \
                                                    sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = INSIDE;             \
                                                    break;                                                                            \
                                            }                                                                                         \
                                            break;                                                                                    \
                                        default: /* if the new value is not zero, then it is a boundary */                            \
                                            sd->shape_next[(d * sd->grid_l[1] + r) * sd->grid_l[0] + c] = BOUNDARY;                   \
                                            break;                                                                                    \
                                    }

#endif // SHAPE_INCL