#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "multigrid.h"
#include <mpi.h>

#define PI 3.141592653589793

int main(int argc, char **argv) 
{
	int nlevel = 15;
	double length_x = 1.0, length_y = 1.0;
    int maxiteration = 20;
	int ngrid;
	int numgrid_x[nlevel], numgrid_y[nlevel]; //Dimension of grid in each level
	int matrixDOF[nlevel]; //Number of elements in each level
	int count, ii, jj, kk, i, iter;  
	double gridsize[nlevel]; //Physical Size of grid
	double *rhs[nlevel], *solution[nlevel]; //Explicit Values of grid
	double coord_x, coord_y;
	clock_t time_start, time_end;
	double time_elapsed;

	char myfilename[256];
	FILE *fp;

	time_start=clock();
    // Beginning of the code
    //
    int wsize, wrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

	ngrid = pow(2,nlevel-1);

    int plevel = 1;

    while(1)
    {
        if (pow(2,plevel) == wsize)
            break;
        plevel += 1;
    }

	numgrid_x[nlevel-1] = ngrid/wsize;
	numgrid_y[nlevel-1] = ngrid;
	gridsize[nlevel-1] = length_y/(double)ngrid;
	matrixDOF[nlevel-1] = numgrid_x[nlevel-1]*numgrid_y[nlevel-1];

    for(i=nlevel-1; i>=plevel+1; i--)
    {
		numgrid_x[i-1] = numgrid_x[i]/2;
		numgrid_y[i-1] = numgrid_y[i]/2;
		gridsize[i-1] = length_y/numgrid_y[i-1];
		matrixDOF[i-1] = numgrid_x[i-1]*numgrid_y[i-1];
    }

    numgrid_x[plevel] *= wsize;
	for(i=plevel; i>=1; i--)
    {
        numgrid_x[i-1] = numgrid_x[i]/2;
        numgrid_y[i-1] = numgrid_y[i]/2;
		gridsize[i-1] = length_y/numgrid_y[i-1];
		matrixDOF[i-1] = numgrid_x[i-1]*numgrid_y[i-1];
    }
    numgrid_x[plevel] /= wsize;

    if (wrank == 0)
        printf("[Multigrid] Geometry and matrix size initialized.\n");

	for(i=0; i<nlevel; i++) 
	{
        if ((i<plevel)&&(wrank==0))
        {
            rhs[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
            solution[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
        }
        else if ((i==plevel)&&(wrank==0))
        {
            rhs[i] = (double*)malloc(wsize * matrixDOF[i] * sizeof(double));
            solution[i] = (double*)malloc(wsize * matrixDOF[i] * sizeof(double));
        }
        else if ((i==(nlevel-1)-5)&&(wrank==0))
        {
            rhs[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
            solution[i] = (double*)malloc(wsize * matrixDOF[i]*sizeof(double));
        }
        else if ((i>plevel)&&(wrank==0))
        {
            rhs[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
            solution[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
        }
        else if (i>=plevel)
        {
            rhs[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
            solution[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
        }
	}

	count = 0;

	for(ii=0; ii<numgrid_x[nlevel-1]; ii++)
	{
		coord_x = (ii+wrank*numgrid_x[nlevel-1])*gridsize[nlevel-1]+0.5*gridsize[nlevel-1];
		for(jj=0; jj<numgrid_y[nlevel-1]; jj++)
		{
			coord_y = jj*gridsize[nlevel-1]+0.5*gridsize[nlevel-1];
			rhs[nlevel-1][count] = sin(coord_x/length_x*PI)* \
					sin(coord_y/length_y*PI);
			count++;
		} 
	}
    if (wrank == 0) {
        printf("[Multigrid] Geometry and rhs constructed.\n");
        printf("[Multigrid] Start solving equations.\n");
    }

	for(iter=0;iter<maxiteration;iter++)
	{
        if (wrank == 0)
            printf("[Multigrid] iteration %d \n",iter+1);

		for(i=nlevel-1; i>plevel; i--) 
		{
            restriction(rhs[i], rhs[i-1], numgrid_x[i-1], numgrid_y[i-1]);
            if (wrank == 0)
            {
                printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", i, i-1);
            }
		}
        // Gather info to Root
        //
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(rhs[plevel], matrixDOF[plevel], MPI_DOUBLE, \
                    rhs[plevel], matrixDOF[plevel], MPI_DOUBLE, 0, \
                    MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (wrank == 0)
        {
            //printf("%p\t%p\n", rhs[plevel], solution[plevel]);
            // Restriction on Root
            for(i=plevel; i>0; i--)
            {
                restriction(rhs[i], rhs[i-1], numgrid_x[i-1], numgrid_y[i-1]);
                printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", i, i-1);
            }

            solution[0][0] = -rhs[0][0]/4.0;
            printf("[Multigrid] Solution at the coarsest level = %f\n",solution[0][0]);

            // Interpolation on Root
            for(i=1; i<plevel; i++)
            {
                interpolation(solution[i-1],solution[i], numgrid_x[i], numgrid_y[i], wrank, wsize);
                printf("[Multigrid] Interpolation of solution done (level %d -> level %d).\n", i-1, i);
	        }
            i = plevel;
            interpolation(solution[i-1],solution[i],numgrid_x[i] * wsize,numgrid_y[i], wrank, wsize);
            printf("[Multigrid] Interpolation of solution done (level %d -> level %d).\n", i-1, i);
        }
        
        // Distribute to the Process
        //
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatter(solution[plevel], matrixDOF[plevel], MPI_DOUBLE, \
                solution[plevel], matrixDOF[plevel], MPI_DOUBLE, \
                0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

	    // Interpolation in Parallel	
        for(i=plevel+1; i<nlevel; i++)
		{
			interpolation(solution[i-1],solution[i], numgrid_x[i], numgrid_y[i], wrank, wsize);
            if (wrank == 0)
                printf("[Multigrid] Interpolation of solution done (level %d -> level %d).\n", i-1, i);
		}

        MPI_Barrier(MPI_COMM_WORLD);
		for(ii=0;ii<matrixDOF[nlevel-1];ii++)
		{
			rhs[nlevel-1][ii] = PI*PI*PI*PI*solution[nlevel-1][ii];
		}
        MPI_Barrier(MPI_COMM_WORLD);
	}

    i=(nlevel-1)-5;
    if (i>=plevel)
    {
        numgrid_x[i] *= wsize;
        MPI_Gather(solution[i], matrixDOF[i], MPI_DOUBLE, \
                    solution[i], matrixDOF[i], MPI_DOUBLE, 0, \
                    MPI_COMM_WORLD);
    }
    if (wrank==0)
    {
        printf("[Multigrid] Final solution obtained.\n");

        sprintf(myfilename,"./solution.dat");
        fp = fopen(myfilename, "wt");

        for(ii=-8; ii<8; ii++)
        {
            coord_x = (ii+numgrid_x[i]/2)*gridsize[i]+0.5*gridsize[i];
            for(jj=-8; jj < 8; jj++)
            {
                coord_y = (jj+numgrid_y[i]/2)*gridsize[i]+0.5*gridsize[i];
                count = (ii+numgrid_x[i]/2)*numgrid_y[i] + jj+numgrid_y[i]/2;
                fprintf(fp,"%12.6f %12.6f %12.6f\n", coord_x, coord_y, solution[i][count]);
            }
        }

        fclose(fp);

        //printf("%p\t%p\n", rhs[plevel], solution[plevel]);
        printf("[Multigrid] Final solution printed.\n");
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

	for(i = 0; i < nlevel; i++) 
	{
        if(wrank==0)
        { 
            free(solution[i]);
            free(rhs[i]);
        }
        else if (i>=plevel)
        {
            free(rhs[i]);
            free(solution[i]);
        }
	}

    if (wrank==0)
        printf("[Multigrid] Memory deallocated.\n");
    MPI_Finalize();
    // End of the Code
	time_end=clock();
	time_elapsed = (double)(time_end - time_start)/CLOCKS_PER_SEC;
    if (wrank==0)
        printf("[Multigrid] Mission completed in %f (secs).\n", time_elapsed);
       
	return 0;
}

