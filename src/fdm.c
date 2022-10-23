#include "fdm.h"
#include "io.h"

#include <mpi.h>

#include <stdlib.h>
#include <math.h>
#include <memory.h>

void DiscretizeSubdomain(Subdomain *sd)
{
    if(sd->grid_g[2] == 0)
    {
        sd->grid_g[0] = sd->dims_l[0] / sd->dx;
        sd->grid_g[1] = sd->dims_l[1] / sd->dy;
        sd->grid_g[2] = 0;
    }
    else
    {
        sd->grid_g[0] = sd->dims_l[0] / sd->dx;
        sd->grid_g[1] = sd->dims_l[1] / sd->dy;
        sd->grid_g[2] = sd->dims_l[2] / sd->dz;
    }
} // end void DiscretizeSubdomain(Subdomain sd)

void DetermineStepSizes(Subdomain *sd)
{
    // step sizes
    sd->dx = (float)sd->dims_g[0] / (float)sd->grid_g[0];
    sd->dy = (float)sd->dims_g[1] / (float)sd->grid_g[1];
    if(sd->n_dims == 3)
    {
        sd->dz = (float)sd->dims_g[2] / (float)sd->grid_g[2];
        sd->mu_z = sd->dt / pow(sd->dz, 2.0);
        if (sd->mu_z > 0.5) { fprintf(stderr, "Step sizes are too large. dz = %f, mu_z = %f > 0.5\n", sd->dx, sd->mu_z); exit(1); }
    }

    // stability factors
    sd->mu_x = sd->dt / pow(sd->dx, 2.0);
    sd->mu_y = sd->dt / pow(sd->dy, 2.0);


    // make sure step sizes are appropriate for stable results
    if (sd->mu_x > 0.5) { fprintf(stderr, "Step sizes are too large. dx = %f, mu_x = %f > 0.5\n", sd->dx, sd->mu_x); exit(1); }
    if (sd->mu_y > 0.5) { fprintf(stderr, "Step sizes are too large. dy = %f, mu_y = %f > 0.5\n", sd->dx, sd->mu_y); exit(1); }
} // end void DetermineStepSizes(Subdomain sd)


