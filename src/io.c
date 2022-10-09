#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ReadInput()
{

}

void WriteData(Subdomain sd)
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
