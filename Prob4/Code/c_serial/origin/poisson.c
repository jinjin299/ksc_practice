#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "matrixconstructor.h"
#include "cgsolver.h"
#include "multigrid.h"

#define PI 3.141592653589793


int main(int argc, char **argv) 
{
    int nlevel, initial_ngrid;
    int *numgrid, *matrixDOF;
    int count, ii, jj, kk, maxiteration;  
    double length, tolerance, *gridsize;
    double **poissonmatrix, **rhs, **solution, **coordinate;
    clock_t time_start, time_end;
    double elapsed_time;

    char myfilename[256];
    FILE *fp;

    nlevel = 3;
    initial_ngrid = 128;
    maxiteration = 20000;
    tolerance = 1.0e-10;
    length = 1.0;
    length = 1.0;

    time_start = clock();

    numgrid = (int*)malloc(nlevel*sizeof(int));
    matrixDOF = (int*)malloc(nlevel*sizeof(int));
    gridsize = (double*)malloc(nlevel*sizeof(double));

    numgrid[nlevel-1] = initial_ngrid;
    gridsize[nlevel-1] = length/initial_ngrid;
    matrixDOF[nlevel-1] = numgrid[nlevel-1]*numgrid[nlevel-1];

    for(ii=nlevel-1; ii>=1; ii--)
    {
	numgrid[ii-1] = numgrid[ii]/2;
	gridsize[ii-1] = length/numgrid[ii-1];
	matrixDOF[ii-1] = numgrid[ii-1]*numgrid[ii-1];
    }

    printf("[Poisson] Geometry and matrix size initialized.\n");

    poissonmatrix = (double**)malloc(nlevel*sizeof(double*));
    rhs = (double**)malloc(nlevel*sizeof(double*));
    solution = (double**)malloc(nlevel*sizeof(double*));
    coordinate = (double**)malloc(nlevel*sizeof(double*));

    for(ii=0; ii<nlevel; ii++) 
    {
	poissonmatrix[ii] = (double*)malloc(matrixDOF[ii]*matrixDOF[ii]*sizeof(double));
	rhs[ii] = (double*)malloc(matrixDOF[ii]*sizeof(double));
	solution[ii] = (double*)malloc(matrixDOF[ii]*sizeof(double));
	coordinate[ii] = (double*)malloc(matrixDOF[ii]*2*sizeof(double));
    }

    for(ii=0; ii<nlevel; ii++)
    {
	count = 0;

	for(jj=1; jj<=numgrid[ii]; jj++)
	{
	    for(kk=1; kk<=numgrid[ii]; kk++)
	    {
		coordinate[ii][2*count] = jj*gridsize[ii]-0.5*gridsize[ii];
		coordinate[ii][2*count+1] = kk*gridsize[ii]-0.5*gridsize[ii];
		count++;
	    } 
	}
    }

    count = 0;

    for(ii=1; ii<=numgrid[nlevel-1]; ii++)
    {
	for(jj=1; jj<=numgrid[nlevel-1]; jj++)
	{
	    rhs[nlevel-1][count] = sin(coordinate[nlevel-1][2*count]/length*PI)* \
				   sin(coordinate[nlevel-1][2*count+1]/length*PI)*gridsize[nlevel-1]*gridsize[nlevel-1];
	    count++;
	} 
    }

    printf("[Poisson] Geometry and rhs constructed.\n");

    for(ii=0; ii<nlevel; ii++) 
	for(jj=0; jj<matrixDOF[ii]*matrixDOF[ii]; jj++)	
	    poissonmatrix[ii][jj] = 0.0;

    for(ii=0; ii<nlevel; ii++)
	construct_poissonmatrix(1, 1, numgrid[ii], numgrid[ii], numgrid[ii], poissonmatrix[ii]);

//        for(jj=0; jj<matrixDOF[0]*matrixDOF[0]; jj++)
//            printf("poisson=%f\n", poissonmatrix[0][jj]);

    printf("[Poisson] Poisson matrix constructed.\n");

    for(ii=0; ii<matrixDOF[nlevel-1]; ii++)
	solution[nlevel-1][ii] = 1.0;

    printf("[Poisson] Start solving equations.\n");

    for(ii=nlevel-1; ii>0; ii--) 
    {
	restriction(rhs[ii], rhs[ii-1], numgrid[ii-1]);

//        for(int aa=0; aa<numgrid[ii-1]*numgrid[ii-1]; aa++)
//                printf("%f\n", rhs[ii-1][aa]);

	printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", ii, ii-1);

	for(jj=0; jj<matrixDOF[ii-1]; jj++)
	    solution[ii-1][jj] = 1.0;
    }

    cgsolver(matrixDOF[0], poissonmatrix[0], rhs[0], solution[0], maxiteration, tolerance);

    for(ii=1; ii<nlevel; ii++)
    {
	interpolation(solution[ii-1], solution[ii], numgrid[ii]);
	printf("[Multigrid] Interpolation of solution done (level %d -> level %d).\n", ii-1, ii);
	cgsolver(matrixDOF[ii], poissonmatrix[ii], rhs[ii], solution[ii], maxiteration, tolerance);
	printf("[Poisson] Solution in level %d obtained.\n", ii);
    }

    printf("[Poisson] Final solution obtained.\n");

    system("rm -rf result\n");
    system("mkdir result\n");

    for(ii=0; ii<nlevel; ii++) 
    {
	sprintf(myfilename,"./result/solution_level%d.dat", ii);
	fp = fopen(myfilename, "wt");

	count = 0;
	for(jj=1; jj<=numgrid[ii]; jj++)
	{
	    for(kk=1; kk<=numgrid[ii]; kk++)
	    {
		fprintf(fp,"%9.5f%10.5f%13.6f\n", coordinate[ii][2*count], coordinate[ii][2*count+1], solution[ii][count]);
		count++;
	    }
	}

	fclose(fp);
    }

    printf("[Poisson] Final solution printed.\n");

    for(ii=0; ii<nlevel; ii++) 
    {
	free(coordinate[ii]);
	free(solution[ii]);
	free(rhs[ii]);
	free(poissonmatrix[ii]);
    }

    free(poissonmatrix);
    free(rhs);
    free(solution);
    free(coordinate);
    free(matrixDOF);
    free(numgrid);
    free(gridsize);

    printf("[Poisson] Memory deallocated.\n");

    time_end = clock();
    elapsed_time = (double)(time_end - time_start)/CLOCKS_PER_SEC;

    printf("[Poisson] Mission completed in %f (secs).\n", elapsed_time);

    return 0;
}


void cgsolver(int size, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance)
{
    int ii, jj, kk;
    double alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0;
    double *p, *Ax, *Ap, *res;

    res = (double*) malloc(size*sizeof(double));
    p   = (double*) malloc(size*sizeof(double));
    Ax  = (double*) malloc(size*sizeof(double));
    Ap  = (double*) malloc(size*sizeof(double));

    multiply(size, matrix, solution, Ax);

    for (ii=0; ii<size; ii++)
    {
	res[ii] = rhs[ii]-Ax[ii];
	p[ii]   = res[ii];

//printf("In CG, %f %f\n", rhs[ii], res[ii]);
    }

    res0tol = innerproduct(res, res, size);

//printf("res0tol=%f\n", res0tol);
//exit(1);

    printf("[CG] Conjugate gradient is started.\n");

    for (ii=0; ii<maxiteration; ii++)
    {
	if ((ii%20==0)&&(ii!=0))
	    printf("[CG] mse %e with a tolerance criteria of %e at %5d iterations.\n", sqrt(temp2/res0tol), tolerance, ii);

	temp1 = innerproduct(res, res, size);
	multiply(size, matrix, p, Ap);
	temp2 = innerproduct(Ap, p, size);

	alpha=temp1/temp2;

	for (jj=0; jj<size; jj++)
	{
	    solution[jj] = solution[jj] + alpha*p[jj];
	    res[jj] = res[jj] - alpha*Ap[jj];
	}

	temp2 = innerproduct(res, res, size);

	if (sqrt(temp2/res0tol) < tolerance)
	    break;

	beta = temp2/temp1;

	for (jj=0; jj<size; jj++)
	    p[jj]= res[jj] + beta*p[jj];

    }

    printf("[CG] Finished with total iteration = %d, mse = %e.\n", (ii+1), sqrt(temp2/res0tol));

    free(res);
    free(p);
    free(Ax);
    free(Ap);
}

double innerproduct(double *x, double *y, int size)
{
    int ii;
    double result;

    result = 0.0;

    for(ii=0; ii<size; ii++)
	result += x[ii]*y[ii];

    return result;
}

void multiply(int size, double *matrix, double *x, double *y)
{
    int ii, jj;

    for (ii=0; ii<size; ii++)       // initialize y
	y[ii]=0.0;

    for (ii=0; ii<size; ii++)
	for (jj=0; jj<size; jj++)
	    y[ii] += matrix[ii*size+jj] * x[jj];
}

