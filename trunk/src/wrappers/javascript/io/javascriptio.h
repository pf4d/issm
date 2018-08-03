/*\file matlabio.h
 *\brief: I/O for ISSM in matlab mode
 */

#ifndef _JAVASCRIPT_IO_H_
#define _JAVASCRIPT_IO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 

#include "../include/javascriptincludes.h"
#include "../../../c/bamg/bamgobjects.h"
#include "../../../c/classes/classes.h"
#include "../../../c/toolkits/toolkits.h"
#include "../../../c/shared/shared.h"

void WriteData(IssmPDouble** pmatrix,int* pnel, int* matrix, int M,int N);
void WriteData(IssmPDouble** pmatrix,int* pM, int* pN, int* matrix, int M, int N);
void WriteData(IssmPDouble** pmatrix,int* pM, int* pN, IssmPDouble* matrix, int M, int N);
void WriteData(IssmPDouble** px, int* pnods, int* vector, int M);
void WriteData(IssmPDouble** px, int* pnods, double* vector, int M);
void WriteData(char** pstring, char* stringin);
void WriteData(IssmPDouble* pdouble, IssmPDouble doublein);
void WriteData(IssmPDouble** pdouble, void*);

void FetchData(char** pstring, char* stringin);
void FetchData(double* pscalar,double scalar);
void FetchData(int* pinteger,int integer);
void FetchData(double** pvector, double* vectorin, int nods);
void FetchData(double** pvector, int* pnods, double* vectorin, int nods);
void FetchData(double **pmatrix, int* pM, int* pN, int* matrixin, int M, int N);
void FetchData(double **pmatrix, int* pM, int* pN, double* matrixin, int M, int N);
void FetchData(int **pmatrix, int* pM, int* pN, int* matrixin, int M, int N);
void FetchData(Contours** pcontours,double* x, double* y, int nods);
void FetchData(Options** poptions,int NRHS, int nrhs, const char* optionname, double optionvalue);
void FetchData(int* pinteger,int integer);

/*Print*/
void ApiPrintf(const char* string);
#endif	/* _IO_H_ */
