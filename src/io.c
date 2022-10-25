#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ReadInput()
{

}

void WriteData(Subdomain sd, int type)
{
    if (type == FDM)
    {
        FILE *fp = fopen("data/data.csv", "a");
        if (fp == NULL) { printf("file not found\n"); exit(1); }

        //fprintf(fp, "%f ", time);
        for (int i = 0; i < sd.grid_g[0] * sd.grid_g[1] * sd.grid_g[2]; i++)
        {
            fprintf(fp, "%f,", sd.u_global[i]);
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
