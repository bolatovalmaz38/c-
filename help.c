#include<stdio.h>
#include<stdlib.h>
#include"mpi.h"
#include"help.h"
#include<math.h>
double f(int i,int j){
	//return 1.0/(i + j + 1);
	return fabs(i-j);
}	
void init_data(double *A,double *B,int n,int my_rank,int p,int max_rows, int flag, char *argv[]){
	int i,j;
	if (flag == 2){
	    	for (i = 0; i < max_rows; i++){
			B[i] = 0.;
			for (j = 0; j < n; j++)
				A[i*n+j] = f(i*p+my_rank,j);
			for (j = 0; j < n; j += 2)
				B[i]+=A[i*n+j];
		}
	}
	else if (flag == 3){
		FILE* file;
		file = fopen(argv[2], "r");
		if (file == NULL)
			printf("Error: file not opened \n");
		double *a,*b;
		a = (double*)malloc(n*n*sizeof(double));
		b = (double*)malloc(n*sizeof(double));
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if((fscanf(file, "%lf", &a[i * n + j])) != 1){
					printf("Error: A[i][j] is incorrect \n");
				}
			}
			if((fscanf(file, "%lf", &b[i])) != 1){
				printf("Error: b[i] is incorrect \n");
			}
		}
		for (i = 0; i < max_rows; i++)
		{
			B[i] = 0.;
			for (j = 0; j < n; j++)
				A[i*n+j] = a[(i*p+my_rank)*n+j];
			for (j = 0; j < n; j += 2)
				B[i] += b[i*p+my_rank];
		}
		free(a);
		free(b);
		fclose(file);
	}
}

void print_data(double *A,double *B,double *v2,int n,int my_rank,int p,int mi,int mj){
	int i,j;
	MPI_Status stat;
	if(!my_rank){
		printf("/*---------------------------------*/\n");
		for(i = 0; i < mi; i++){
			if(i%p == 0){
				for(j = 0; j < mj; j++)
					printf("[%3.4lf]", A[(i/p)*n+j]);
				printf("|[%3.4lf]\n", B[i]);
			}
			else{
				MPI_Recv(v2, mj + 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &stat);
				for(j = 0; j < mj; j++)
					printf("[%3.4lf]",v2[j]);
				printf("|[%3.4lf]\n",v2[mj]);
			}
		}
		printf("/*---------------------------------*/\n");	
	}
	else{
		for(i = 0; i * p + my_rank < mi; i++){
			for(j = 0; j < mj; j++)
				v2[j] = A[i*n+j];
			v2[mj] = B[i];
			MPI_Send(v2, mj + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}			
