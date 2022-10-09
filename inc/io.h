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
void WriteData(Subdomain subdomain, int type); 

#endif // IO_INCL