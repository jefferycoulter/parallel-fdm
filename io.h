#ifndef IO_INCL
#define IO_INCL

#include <stdio.h>

/**
 * @brief read user input and assign to variables
 */
void ReadInput();

/**
 * @brief write each entire iteration result to a single line
 * @param data the data to be written
 */
void WriteData(float *data, int cells[3]); 

/**
 * @brief write data in matrix order, i.e. first row is first row of matrix, 
 * second row is second row of matrix, and so on
 * @param data the data to be written
 */
void WriteDataOrdered(float *data, float time, int cells[3]); 

#endif // IO_INCL