#include "fdm.h"
#include "io.h"

#include <stdlib.h>
#include <math.h>
#include <memory.h>

void ComputeFD(Data data, StepSize dX, int time_steps, int cells[3])
{
    WriteData(data.now, cells); // save initial conditions
    for (int t = 1; t < time_steps+1; t ++)
    {
        FTCS(data, dX, cells);
        WriteData(data.next, cells); // save each iteration
    }
}

void FTCS(Data data, StepSize dX, int cells[3])
{
    float mu_x = dX.dt / pow(dX.dx, 2.0);
    float mu_y = dX.dt / pow(dX.dy, 2.0);
    // possible that only 2D is considered, so check
    float mu_z = dX.dz == 0.0 ?  0.0 : dX.dt / pow(dX.dz, 2.0);

    // compute finite difference
    for (int x = 1; x < cells[0] - 1; x++)
    {
        for (int y = 1; y < cells[1] - 1; y++)
        {
            if (cells[2] == 0.0) // 2D problem
            {
                data.next[x * cells[1] + y] = (1 - 2 * mu_x - 2 * mu_y) * data.now[x * cells[1] + y]          \
                    + mu_x * (data.now[x * cells[1] + y + cells[0]] + data.now[x * cells[1] + y - cells[0]])  \
                    + mu_y * (data.now[x * cells[1] + y + 1] + data.now[x * cells[1] + y - 1]);
            }
            else
            {
                for (int z = 1; z < cells[2] - 1; z++)
                {
                    data.next[(x * cells[1] + y) * cells[2] + z] = (1 - 2 * mu_x - 2 * mu_y - 2 * mu_z) * data.now[(x * cells[1] + y) * cells[2] + z]                   \
                    + mu_x * (data.now[(x * cells[1] + y + cells[0]) * cells[2] + z] + data.now[(x * cells[1] + y - cells[0]) * cells[2] + z])    \
                    + mu_y * (data.now[(x * cells[1] + y + 1) * cells[2] + z] + data.now[(x * cells[1] + y - 1) * cells[2] + z])                  \
                    + mu_z * (data.now[(x * cells[1] + y) * cells[2] + z + cells[0] * cells[1]] + data.now[(x * cells[1] + y) * cells[2] + z - cells[0] * cells[1]]);
                }
            }
        }
    }

    // copy new data to old data array for next iteration
    memcpy(data.now, data.next, cells[0] * cells[1] * sizeof((*data.now)));
}

void SetBoundaryConditions(Data data, int cells[3])
{
    if (cells[2] == 0.0) // 2D problem
    {
        for (int x = 0; x < cells[0]; x++)
        {
            data.now[x * cells[1]] = 10; // left side
            data.next[x * cells[1]] = 10;
            if (x == 0) // fill in top
            {
                for (int y = 1; y < cells[1]; y++)
                {
                    data.now[x * cells[1] + y] = 10;
                    data.next[x * cells[1] + y] = 10;
                }
            }
            else if (x == cells[0] - 1) // fill in bottom
            {
                for (int y = 1; y < cells[1]; y++)
                {
                    data.now[x * cells[1] + y] = 10;
                    data.next[x * cells[1] + y] = 10;
                }
            }
            else // fill in the rest of the terms on the right side
            {
                data.now[x*cells[1] + cells[1] - 1] = 10;
                data.next[x*cells[1] + cells[1] - 1] = 10;
            }
        }
    }
    else // 3D problem
    {
        for (int x = 0; x < cells[0]; x++)
        {
            data.now[x * cells[1]] = 10; // left side bottom edge
            data.next[x * cells[1]] = 10;
            if (x == 0) // fill in top
            {
                for (int y = 1; y < cells[1]; y++)
                {
                    data.now[x * cells[1] + y] = 10;
                    data.next[x * cells[1] + y] = 10;
                }
            }
            else if (x == cells[0] - 1) // fill in bottom
            {
                for (int y = 1; y < cells[1]; y++)
                {
                    data.now[x * cells[1] + y] = 10;
                    data.next[x * cells[1] + y] = 10;
                }
            }
            else // fill in the rest of the terms on the right side
            {
                data.now[x*cells[1] + cells[1] - 1] = 10;
                data.next[x*cells[1] + cells[1] - 1] = 10;
            }
        }
    }

}

void SetInitialConditions(float *data, int cells[3])
{
    for (int x = 1; x < cells[0] - 1; x++)
    {
        for (int y = 1; y < cells[1] - 1; y++)
        {
            if (cells[2] == 0) // 2D problem
            {
                data[x * cells[1] + y] = -10;
            }
            else // 3D problem
            {
                for (int z = 1; z < cells[2] - 1; z++)
                {
                    data[(x * cells[1] + y) * cells[2] + z] = -10;
                }
            }
        }
    }
}

void DiscretizeSubdomain(Dimension dim, StepSize dX, int cells[3])
{
    if(cells[2] == 0)
    {
        cells[0] = dim.x / dX.dx;
        cells[1] = dim.y / dX.dy;
    }
    else
    {
        cells[0] = dim.x / dX.dx;
        cells[1] = dim.y / dX.dy;
        cells[2] = dim.z / dX.dz;
    }
}

void GenerateSpatialArrays(Data *data, int cells[3])
{
    if (cells[2] == 0.0) // 2D problem
    {
        data->next = (float*)malloc(cells[0] * cells[1] * sizeof(float));
        data->now = (float*)malloc(cells[0] * cells[1] * sizeof(float));
        memset(data->next, 0, cells[0] * cells[1] * sizeof((*data->next)));
        memset(data->now, 0, cells[0] * cells[1] * sizeof((*data->now)));
    }
    else // 3D problem
    {
        data->next = (float*)malloc(cells[0] * cells[1] * cells[2] * sizeof(float));
        data->now = (float*)malloc(cells[0] * cells[1] * cells[2] * sizeof(float));
        memset(data->next, 0, cells[0] * cells[1] * cells[2] * sizeof((*data->next)));
        memset(data->now, 0, cells[0] * cells[1] * cells[2] * sizeof((*data->now)));
    }
}



