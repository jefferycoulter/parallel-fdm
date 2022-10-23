#ifndef IO_INCL
#define IO_INCL

#include "parallel.h"
#include <stdio.h>

enum DataType{FDM, Shape};

/**
 * @brief read user input and assign to variables
 */
void ReadInput();

/**
 * @brief write each entire iteration result to a single line
 * @param Subdomain the data to be written
 */
#define WriteData(sd, data, type)                                                       \
    if (type == FDM)                                                                    \
    {                                                                                   \
        FILE *fp = fopen("data/data.csv", "a");                                         \
        if (fp == NULL) { printf("file not found\n"); exit(1); }                        \
        for (int i = 0; i < sd->grid_g[0] * sd->grid_g[1] * sd->grid_g[2]; i++)         \
        {                                                                               \
            fprintf(fp, "%p,", (void*)&data);                                           \
        }                                                                               \
        fprintf(fp, "\n");                                                              \
        fclose(fp);                                                                     \
    }                                                                                   \
    else if (type == Shape)                                                             \
    {                                                                                   \
        FILE *fp = fopen("data/shape.csv", "a");                                        \
        if (fp == NULL) { printf("file not found\n"); exit(1); }                        \
        for (int i = 0; i < sd->grid_g[0]; i++)                                         \
        {                                                                               \
            for (int j = 0; j < sd->grid_g[1]; j++)                                     \
            {                                                                           \
                for (int k = 0; k < sd->grid_g[2]; k++)                                 \
                {                                                                       \
                    fprintf(fp, "%d,", data[i][j][k]);                                  \
                }                                                                       \
            }                                                                           \
                                                                                        \
        }                                                                               \
        fprintf(fp, "\n");                                                              \
        fclose(fp);                                                                     \
    }                                                                                   \
    else                                                                                \
    {                                                                                   \
        fprintf(stderr, "Choose valid type of data to write: FDM or Shape\n");          \
        exit(1);                                                                        \
    }
//void WriteData(Subdomain subdomain, int type); 

#endif // IO_INCL