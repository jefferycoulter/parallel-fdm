#ifndef SHAPE_INCL
#define SHAPE_INCL

#include "fdm.h"
#include "parallel.h"

/**
 * @brief create a shape array.  in the case of irregular geometry, this corresponds to which locations
 * are inside the area of interest (i.e. where FDM computation occurs), which locations are boundaries,
 * and which locations are outside the area of interest
 * @param subdomain subdomain to create shape array in
 * @param radius radius of domain of interest
 */
// void CreateShapeArray(Subdomain *subdomain, float radius);

/**
 * @brief apply a laplace filter (discrete laplace operator) to detect the boundary of the area
 * of interest (i.e. where FDM computation occurs)
 * @param subdomain subdomain to apply filter on
 */
// void ApplyLaplaceFilter(Subdomain *subdomain);

/**
 * @brief shift the coordinates in a given subdomain so that the center of the global subdomain (dim_X / 2, dim_y / 2) 
 * is moved to the origin (0,0).  this is used for determining boundaries in the case of irregular geometry
 * @param subdomain subdomain to shift
 * @param radius the radius of the irregular geometry
 */
#define CoordShift(sd, s0, r)                                                                                                           \
    for (int i = sd->bounds_l[0]; i < sd->bounds_l[1]; i ++)                                                                            \
    {                                                                                                                                   \
        for (int j = sd->bounds_l[2]; j < sd->bounds_l[3]; j++)                                                                         \
        {                                                                                                                               \
            for (int k = sd->bounds_l[4]; k < sd->bounds_l[5]; k++)                                                                     \
            {                                                                                                                           \
                if (pow(i - (sd->grid_g[0] / 2), 2.0) + pow(j - (sd->grid_g[1] / 2), 2.0) + pow(k - (sd->grid_g[2] / 2), 2.0) <= r*r)   \
                {                                                                                                                       \
                    s0[i - sd->bounds_l[0]][j - sd->bounds_l[2]][k - sd->bounds_l[4]] = 1;                                              \
                }                                                                                                                       \
                else                                                                                                                    \
                {                                                                                                                       \
                    s0[i - sd->bounds_l[0]][j - sd->bounds_l[2]][k - sd->bounds_l[4]] = 0;                                              \
                }                                                                                                                       \
            }                                                                                                                           \
        }                                                                                                                               \
    }
//void CoordShift(Subdomain *subdomain, float radius);

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
#define AssignValue(s, s0, i, j, k)    switch (s[i][j][k])                                                                            \
                                    {                                                                                                 \
                                        case 0: /* if the new value is zero, then the previous value was either INSIDE or OUTSIDE*/   \
                                            switch (s0[i][j][k])                                                                      \
                                            {                                                                                         \
                                                case Outside: /* if previous value was OUTSIDE, then new value is still outside */    \
                                                    s[i][j][k] = Outside;                                                             \
                                                    break;                                                                            \
                                                                                                                                      \
                                                case Inside: /* if previous value was INSIDE, then switch it back to inside */        \
                                                    s[i][j][k] = Inside;                                                              \
                                                    break;                                                                            \
                                            }                                                                                         \
                                            break;                                                                                    \
                                        default: /* if the new value is not zero, then it is a boundary */                            \
                                            s[i][j][k] = Boundary;                                                                    \
                                            break;                                                                                    \
                                    }

/**
 * @brief apply a laplace filter (discrete laplace operator) to detect the boundary of the area
 * of interest (i.e. where FDM computation occurs)
 * @param sd subdomain to apply filter on
 * @param s resulting shape array
 * @param s0 this contains the initial shape prior
 * @param r radius of circle
 */
#define ApplyLaplaceFilter(sd, s, s0, mu, r)                \
    for (int i = 0; i < sd->grid_l[0]; i++)                 \
    {                                                       \
        for (int j = 0; j < sd->grid_l[0]; j++)             \
        {                                                   \
            for (int k = 0; k < sd->grid_l[0]; k++)         \
            {                                               \
                InteriorFD(s, s0, mu, i, j, k)              \
                AssignValue(s, s0, i, j, k)                 \
            }                                               \
        }                                                   \
    }
//void ApplyLaplaceFilter(Subdomain *subdomain);

/**
 * @brief create a shape array.  in the case of irregular geometry, this corresponds to which locations
 * are inside the area of interest (i.e. where FDM computation occurs), which locations are boundaries,
 * and which locations are outside the area of interest
 * @param subdomain subdomain to create shape array in
 * @param radius radius of domain of interest
 */
#define CreateShapeArray(sd, s, s0, mu, r)  \
    CoordShift(sd, s0, r)                   \
    ApplyLaplaceFilter(sd, s, s0, mu, r);   \
//void CreateShapeArray(Subdomain *subdomain, float radius);

#endif // SHAPE_INCL