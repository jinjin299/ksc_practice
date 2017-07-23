#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kmeans.h"
#include <mpi.h>
void assignment_step(const PPOINT kmeans, const PPOINT pt, int *kindex, int n_PT)
{
	double d1, d2;

	for(int i=0; i<n_PT; i++) {
		kindex[i] = 0;
		d1 = distance(pt[i], kmeans[0]);
		for(int j=1; j<N_K; j++) {
			d2 = distance(pt[i], kmeans[j]);
			if( d1 > d2 ) {
				d1 = d2;
				kindex[i] = j;
			}
		}
	}
}

void update_step(PPOINT kmeans, const PPOINT pt, const int *kindex, int n_PT)
{
	int i, idx;
	int num_pt[N_K];
	int numr_pt[N_K];
    double kx[N_K], ky[N_K];
    double krx[N_K], kry[N_K];

	for(i=0; i<N_K; i++) {
        kx[i] = ky[i] = 0.0;
		//kmeans[i].x = kmeans[i].y = 0.0;
		num_pt[i] = 0;
	}
	for(i=0; i<n_PT; i++) {
		idx = kindex[i];
        kx[idx] += pt[i].x;
        ky[idx] += pt[i].y;
		//kmeans[idx].x += pt[i].x;
		//kmeans[idx].y += pt[i].y;
		num_pt[idx]++;
	}
    // Reduction of kmeans & num_pt
    MPI_Allreduce(kx, krx, N_K, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(ky, kry, N_K, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(num_pt, numr_pt, N_K, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for(i=0; i<N_K; i++) {
		kmeans[i].x = krx[i]/numr_pt[i];
		kmeans[i].y = kry[i]/numr_pt[i];
	}
}
