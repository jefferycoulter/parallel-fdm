#ifndef FDM_INCL
#define FDM_INCL

#include "io.h"
#include "parallel.h"
#include "shape.h"

#include <mpi.h>

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <math.h>

/**
 * @brief location of fdm grid.  outside region of interest (i.e. no FDM computations here), 
 * region of interest (i.e. where FDM computations occur), or boundary of region of interestion
 */
enum Location{Outside, Inside, Boundary};

/**
 * @brief boundary conditions that can be used for fdm
 */
enum BC{Dirichlet, VonNeumann};

/**
 * @brief stability factors
 */
typedef struct 
{
    float x, y ,z;
} Mu;

/**
 * @brief compute interior points. same for both dirichlet and von neumann
 * @param subdomain the subdomain on which to compute finite differences
 * @param i x index
 * @param j y index
 * @param z index (zero if 2D)
 */
#define InteriorFD(u, u0, mu, i, j, k)                                          \
    if (k == 0)                                                                 \
    {                                                                           \
        u[i][j][k] = (1 - 2 * mu.x - 2 * mu.y) * u0[i][j][k]                    \
                            + mu.x * (u0[i + 1][j][k] + u0[i - 1][j][k])        \
                            + mu.y * (u0[i][j + 1][k] + u0[i][j - 1][k]);       \
    }                                                                           \
    else                                                                        \
    {                                                                           \
            u[i][j][k] = (1 - 2 * mu.x - 2 * mu.y - 2 * mu.z) * u0[i][j][k]     \
                                + mu.x * (u0[i + 1][j][k] + u0[i - 1][j][k])    \
                                + mu.y * (u0[i][j + 1][k] + u0[i][j - 1][k])    \
                                + mu.z * (u0[i][j][k + 1] + u0[i][j][k - 1]);   \
    }
//void InteriorFD(Subdomain *subdomain, int i, int j, int k);

/**
 * @brief compute boundary terms depending on whether dirichlet or von neuman is being used
 * @param subdomain the subdomain on which to compute finite differences
 * @param i x index
 * @param j y index
 * @param z index (zero if 2D)
 */
#define BoundaryFD(u, u0, bc, i, j, k)  \
    switch (bc)                         \
    {                                   \
        case Dirichlet:                 \
        {                               \
            break;                      \
        }                               \
        case VonNeumann:                \
        {                               \
            /* need to implement */     \
            break;                      \
        }                               \
    }
//void BoundaryFD(Subdomain *subdomain, int bc, int i, int j, int k);

/**
 * @brief macro to decide if the point is interior of the domain or a boundary point
 * @param sd subdomain
 * @param bc boundary conditions
 * @param i x index
 * @param j y index
 * @param k z index
 */
#define InteriorOrBoundary(u, u0, s, bc, i, j, k)   \
    switch (s[i][j][k])                             \
    {                                               \
    case Outside:                                   \
        break;                                      \
    case Inside:                                    \
        InteriorFD(u, u0, i, j, k)                  \
        break;                                      \
    case Boundary:                                  \
        BoundaryFD(u, u0, bc, i, j, k)              \
        break;                                      \
    } 

/**
 * @brief compute a single temporal iteration of forward-time central-space (FTCS) finite difference 
 * scheme for a parabolic PDE (2D diffusion equation right now)
 * @param subdomain the subdomain on which to compute finite differences
 * @param bc boundary conditions (Dirichlet or VonNeumann)
 */
#define FTCS(sd, u, u0, bc)                                 \
    for (int i = 0; i < sd->grid_l[0]; i++)                 \
    {                                                       \
        for (int j = 0; j < sd->grid_l[0]; j++)             \
        {                                                   \
            for (int k = 0; k < sd->grid_l[0]; k++)         \
            {                                               \
                InteriorOrBoundary(u, u0, s, bc, i, j, k)   \
            }                                               \
        }                                                   \
    }
//void FTCS(Subdomain *subdomain, int bc);

/**
 * @brief run finite difference method for a given number of iterations
 * @param subdomain the subdomain on which to compute finite differences
 * @param bc boundary conditions (Dirichlet or VonNeumann)
 * @param time_steps max number of iterations for computations
 */
#define ComputeFD(sd, u, u0, bc, ts)                                                    \
    MPI_Gather(u0, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,            \
                   sd->u_g, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,   \
                   ROOT, MPI_COMM_WORLD);                                               \
    WriteData(*sd, FDM); /* save initial conditions */                                  \
    for (int t = 1; t < ts+1; t ++)                                             \
    {                                                                                   \
        FTCS(sd, u, u0, bc); /* compute finite differences */                           \
                                                                                        \
        /* send the local computations to ROOT process */                               \
        MPI_Gather(u, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,         \
                   sd->u_g, sd->grid_l[0] * sd->grid_l[1] * sd->grid_l[2], MPI_FLOAT,   \
                   ROOT, MPI_COMM_WORLD);                                               \
        if (sd->rank == ROOT)                                                           \
        {                                                                               \
            WriteData(*sd, FDM); /* save each iteration */                              \
        }                                                                               \
    } 
//void ComputeFD(Subdomain *subdomain, int bc, int time_steps);

/**
 * @brief apply boundary conditions to the domain
 * @param subdomain the subdomain in which to set boundary conditions
 */
#define SetBoundaryConditions(sd, u, s, bc)             \
    for (int i = 0; i < sd->grid_l[0]; i++)             \
    {                                                   \
        for (int j = 0; j < sd->grid_l[1]; j++)         \
        {                                               \
            for (int k = 0; k < sd-<grid_l[2]; k++)     \
            {                                           \
                switch (s[i][j][k])                     \
                {                                       \
                    case BOUNDARY:                      \
                        u[i][j][k] = bc;                \
                        break;                          \
                    default:                            \
                        break;                          \
                }                                       \
            }                                           \
        }                                               \
    }
//void SetBoundaryConditions(Subdomain *subdomain);

/**
 * @brief apply initial conditions to the domain
 * @param subdomain the subdomain in which to set initial conditions
 */
#define SetInitialConditions(sd, u, s, ic)              \
    for (int i = 0; i < sd->grid_l[0]; i++)             \
    {                                                   \
        for (int j = 0; j < sd->grid_l[1]; j++)         \
        {                                               \
            for (int k = 0; k < sd-<grid_l[2]; k++)     \
            {                                           \
                switch (s[i][j][k])                     \
                {                                       \
                    case BOUNDARY:                      \
                        u[i][j][k] = ic;                \
                        break;                          \
                    default:                            \
                        break;                          \
                }                                       \
            }                                           \
        }                                               \
    }
//void SetInitialConditions(Subdomain *subdomain);

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

void AllocateLocalArray();

void AllocateGlobalArray();

#endif // FDM_INCL