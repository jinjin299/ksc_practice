#ifndef CGSOLVER_H_
#define CGSOLVER_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define Printf if((wrank==0)) printf

void cgsolver(int size, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance, int wrank, int wsize);
double innerproduct(double *x, double *y, int size);
void multiply(int size, double *matrix, double *x, double *y, int sdex, int rsize);

#endif
