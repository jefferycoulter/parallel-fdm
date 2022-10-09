#ifndef IO_INCL
#define IO_INCL

#include "parallel.h"
#include <stdio.h>

/**
 * @brief read user input and assign to variables
 */
void ReadInput();

/**
 * @brief write each entire iteration result to a single line
 * @param data the data to be written
 */
void WriteData(Subdomain subdomain); 

#endif // IO_INCL