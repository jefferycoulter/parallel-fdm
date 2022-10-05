#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ReadInput()
{

}

void WriteData(float *data, int cells[3])
{
    FILE *fp = fopen("data/data.csv", "a");
    if (fp == NULL) { printf("file not found\n"); exit(1); }

    //fprintf(fp, "%f ", time);
    for (int i = 0; i < cells[0] * cells[1]; i++)
    {
        fprintf(fp, "%f,", data[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void WriteDataOrdered(float *data, float time, int cells[3])
{
    FILE *fp = fopen("data/data.txt", "a");
    if (fp == NULL) { printf("file not found\n"); exit(1); }

    fprintf(fp, "%f ", time);
    for (int x = 0; x < cells[0]; x++)
    {
        for (int y = 0; y < cells[1]; y++)
        {
            fprintf(fp, "%f ", data[x * cells[0] + y]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}