#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<sys/types.h>
#include"mpi.h"
#include"sol.h"
#include"help.h"
#include"mspin.h"
#define N_MAX 8
int main(int argc,char* argv[]){
    	int my_rank;
    	int p;
    	int i, n,mi;
    	mspin *spin;
    	int max_rows;
    	double t;
    	double *A;
    	double *B;
    	double *v2, *v, *v3;
    	MPI_Init(&argc,&argv);
    	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	if (argc > 1 && argc < 4)
	{
		sscanf(argv[1], "%d", &n);
		if (n <= 0) {
			printf("n <= 0\n");
			return -2;
		}
	}
	else{
		printf("Incorrect mode!\n");
		return -3;
	}
    	max_rows = 0;
    	for (mi = my_rank; mi < n; mi += p) 
		max_rows++;
    	A = (double*)malloc(n * max_rows * sizeof(double));
    	B = (double*)malloc(max_rows * sizeof(double));
    	v2 = (double*)malloc((n+1) * sizeof(double));
    	v = (double*)malloc(n * sizeof(double));
    	v3 = (double*)malloc(n * sizeof(double));
    	spin = (mspin*)malloc(sizeof(mspin));
    	if(!(A&&B)){
    		printf("Not enough memory!\n");
    	  	MPI_Abort(MPI_COMM_WORLD,1);
    	}
    	mi=(n>N_MAX)?N_MAX:n;
	init_data(A,B,n,my_rank,p,max_rows,argc,argv);
    	MPI_Barrier(MPI_COMM_WORLD);
	for (i = 0; i < n; i++){
		v[i] = 0; 
		v3[i] = 0;
	} 
	if (n > 1){
	t = MPI_Wtime();
    	Solve_it(A,B,v2,n,spin,my_rank,p,max_rows);
	for (i = 0; i < max_rows; i++)
		v[i*p+my_rank] = B[i]; 
	MPI_Allreduce(v, v3, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    	MPI_Barrier(MPI_COMM_WORLD);
    	init_data(A,B,n,my_rank,p,max_rows,argc,argv);
	Norma(A,B,n,v3,my_rank,max_rows);
	t = MPI_Wtime() - t;
	printf("process %d time:  %lf\n", my_rank, t);
	if (my_rank == 0){
		printf("x = [ ");
		for (i = 0; i < mi; i++)
			printf("%lf ", v3[i]);
		puts("]\n");
	}
  	free(A);
   	free(B);
   	free(v2);
   	free(v);
   	free(v3);
	}
	else{
		if (A[0] > 1e-16)//не ноль
			printf ("x = %lf\n", B[0] /= A[0]);
		else
			printf ("x \\in \\mathbb{R}\n");
	}	
    	MPI_Finalize();
    	return 0;
}
