#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ReadInput(Subdomain *sd)
{
    FILE *fp = fopen("config", "r");
    if (fp == NULL) { printf("file not found\n"); exit(1); }
    fscanf(fp, "ndims = %d\n", &(*sd).n_dims);
    fscanf(fp, "dimX = %d\n", &(*sd).dims_g[0]);
    fscanf(fp, "dimY = %d\n", &(*sd).dims_g[1]);
    fscanf(fp, "dimZ = %d\n", &(*sd).dims_g[2]);
    fscanf(fp, "Nx = %d\n", &(*sd).grid_g[0]);
    fscanf(fp, "Ny = %d\n", &(*sd).grid_g[1]);
    fscanf(fp, "Nz = %d\n", &(*sd).grid_g[2]);
    fscanf(fp, "timestep = %f\n", &(*sd).dt);
    fclose(fp);
}

void WriteData(Subdomain sd, int type)
{
    if (type == FDM)
    {
        // species 1
        FILE *fp = fopen("data/data_u.csv", "a");
        if (fp == NULL) { printf("file not found\n"); exit(1); }
        for (int i = 0; i < sd.grid_g[0] * sd.grid_g[1] * sd.grid_g[2]; i++)
        {
            fprintf(fp, "%f,", sd.u_global[i]);
        }
        fprintf(fp, "\n");
        fclose(fp);

        // species 2
        fp = fopen("data/data_v.csv", "a");
        if (fp == NULL) { printf("file not found\n"); exit(1); }
        for (int i = 0; i < sd.grid_g[0] * sd.grid_g[1] * sd.grid_g[2]; i++)
        {
            fprintf(fp, "%f,", sd.v_global[i]);
        }
        fprintf(fp, "\n");
        fclose(fp);

        // species 3
        fp = fopen("data/data_uv.csv", "a");
        if (fp == NULL) { printf("file not found\n"); exit(1); }
        for (int i = 0; i < sd.grid_g[0] * sd.grid_g[1] * sd.grid_g[2]; i++)
        {
            fprintf(fp, "%f,", sd.uv_global[i]);
        }
        fprintf(fp, "\n");
        fclose(fp);
    }
    else if (type == Shape)
    {
        FILE *fp = fopen("data/shape.csv", "a");
        if (fp == NULL) { printf("file not found\n"); exit(1); }
        for (int i = 0; i < sd.grid_g[0] * sd.grid_g[1] * sd.grid_g[2]; i++)
        {
            fprintf(fp, "%d,", sd.shape_g[i]);
        }
        fprintf(fp, "\n");
        fclose(fp);
    }
    else
    {
        fprintf(stderr, "Choose valid type of data to write: FDM or Shape\n");
        exit(1);
    }
}
