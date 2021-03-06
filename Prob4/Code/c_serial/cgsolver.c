#include "cgsolver.h"
#include <mpi.h>
#define Printf if((wrank)==0) printf

void cgsolver(int size, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance, int wsize, int wrank)
{
    int ii, jj, kk;
    double alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0;
    double *p, *Ax, *Ap, *res;

    int rnum=size/wsize;
    res = (double*) malloc(rnum*sizeof(double));
    p   = (double*) malloc(rnum*sizeof(double));

    rest = (double*) malloc(rnum*sizeof(double));
    pt   = (double*) malloc(rnum*sizeof(double));
    Ax  = (double*) malloc(rnum*sizeof(double));
    Ap  = (double*) malloc(rnum*sizeof(double));
    Axt  = (double*) malloc(rnum*sizeof(double));
    Apt  = (double*) malloc(rnum*sizeof(double));
    
    multiply(pnum, size, matrix, solution, Ax);

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
            printf("[CG] mse %e with a tolerance criteria\
                    of %e at %5d iterations.\n", sqrt(temp2/res0tol), tolerance, ii);
        temp1 = innerproduct(res, res, size);
        //Reduction
        //
        multiply(size, matrix, p, Ap);

        // Share and Reduce in each part

        temp2 = innerproduct(Ap, p, size);
        //Reduction

        alpha=temp1/temp2;
        //shared data

        for (jj=0; jj<size; jj++)
        {
            solution[jj] = solution[jj] + alpha*p[jj];
            res[jj] = res[jj] - alpha*Ap[jj];
        }

        temp2 = innerproduct(res, res, size);
        //Reduction

        if (sqrt(temp2/res0tol) < tolerance)
        {
            // Wrap-UP
            break;
        }

        beta = temp2/temp1;
        //independently
        
        for (jj=0; jj<size; jj++)
            p[jj]= res[jj] + beta*p[jj];
        //Indiv

    }
    // Res, p, solution

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

void multiply(int pnum, int size, double *matrix, double *x, double *y)
{
    int ii, jj;

    for (ii=0; ii<size; ii++)       // initialize y
        y[ii]=0.0;

    for (ii=0; ii<size; ii++)
        for (jj=0; jj<pnum; jj++)
            y[ii] += matrix[ii*size+jj] * x[jj];
}

