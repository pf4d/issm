/*!\file:  AverageOntoPartitionx.h
 * \brief header file for averaging  dakota responses onto a vertices partition
 */ 

#ifndef _AVERAGEONTOPARTITIONXX_H
#define _AVERAGEONTOPARTITIONXX_H

#include "../../classes/classes.h"

void AverageOntoPartitionx(double** average, Elements* elements, Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,double* vertex_response);

#endif  /* _AVERAGEONTOPARTITIONXX_H */
