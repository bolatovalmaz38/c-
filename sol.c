#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"mpi.h"
#include"mspin.h"
#include"sol.h"
#include"help.h"
double kor;
void Norma(double *a, double *b, int n, double *x, int r, int row){
	double s, t = 0, g = 0;
	int i,j;
	for (i = 0; i < row; i++){
		s = 0;
		for (j = 0; j < n; j++){
			s += a[i*n+j]*x[j];
		}
		s -= b[i];
		t += s*s;
	}
	MPI_Reduce(&t, &g, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (r == 0)
		printf("Norma = %e\n", sqrt(g));
}
int spin(double *v1, double *v2, int k, mspin *s){
  	kor = sqrt(v1[k] * v1[k] + v2[k] * v2[k]);
  	s->cosf = v1[k] / kor;
  	s->sinf = -v2[k] / kor;
  	s->k=k;
  	s->kor=kor;
  	return 0;
}
void rotate(double *v1, double *v2, double *w1, double *w2, int n, int k, mspin *ms){
	int i;
  	spin(v1, v2, k, ms);
  	for (i = 0; i < ms->k; i++){
    		kor = v1[i];
    		v1[i] = v1[i] * ms->cosf - v2[i] * ms->sinf;
    		v2[i] = kor * ms->sinf + v2[i] * ms->cosf;
  	}
  	v1[ms->k] = ms->kor;
  	v2[ms->k] = 0.;
  	for (i = ms->k + 1; i < n; i++){
    		kor = v1[i];
    		v1[i] = v1[i] * ms->cosf - v2[i] * ms->sinf;
    		v2[i] = kor * ms->sinf + v2[i] * ms->cosf;
  	}
  	kor = w1[0];
  	w1[0] = w1[0] * ms->cosf - w2[0] * ms->sinf;
  	w2[0] = kor * ms->sinf + w2[0] * ms->cosf;
}
void BackGauz(double *A, double *B, double *v2, int n, int my_rank, int p, int max_rows){
	int i,j,frst;
	for (i = 0; i < max_rows; i++){
		B[i] /= A[i * n + i * p + my_rank];
		for (j = n - 1; j >= i * p + my_rank; j--)
			A[i * n + j] /= A[i * n + i * p + my_rank];
	}		
	for (i = n - 1; i >= 1; i--){
		MPI_Barrier(MPI_COMM_WORLD);
		if (my_rank <= i % p)
			frst = i / p;
		else if (i >= p)
			frst = i / p - 1;	
		else{
			MPI_Bcast(v2 + n, 1, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			continue;
		}	
		if(my_rank == i % p){
			MPI_Bcast(B + i / p, 1, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
			for(j = i / p - 1; j >= 0; j--)
				B[j] -= B[i / p] * A[j * n + i];
		}
		else{
			MPI_Bcast(v2 + n, 1, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			for(j = frst; j >= 0; j--)
				B[j] -= v2[n] * A[j*n+i];
		}
	}
}	
void Solve_it(double *A,double *B,double* v2, int n, mspin* ms, int my_rank, int p,int max_rows){
	int i,j,frst, flag;
	MPI_Status status;
	for (i = 0; i < n - 1; i++){
		MPI_Barrier(MPI_COMM_WORLD);
		if (my_rank >= i % p){
			if(i / p >= max_rows)
				continue;
			frst = i / p;
		}	
		else if (i / p + 1 < max_rows)
			frst = i / p + 1;	
		else
			continue;
		for (j = frst + 1; j < max_rows; j++)
			rotate(A + frst * n, A + j * n, B + frst, B + j, n, i, ms);
		if (my_rank == i % p){
			if(i >= n - p + 1){
				flag = (n - i) % p;
				flag = (my_rank + flag) % p;
			}	
			else 
				flag = my_rank;	
			for (j = (my_rank+1) % p; j != flag; j = (j + 1) % p){
	    			MPI_Recv(v2,n+1,MPI_DOUBLE,j,0,MPI_COMM_WORLD,&status);
				rotate(A + frst * n, v2, B + frst, v2 + n, n, i, ms);
				MPI_Send(v2, n + 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
			}		
				
		}
		else{
			for (j = 0; j < n; j++)
				v2[j] = A[frst * n + j];
			v2[n] = B[frst];
			MPI_Send(v2, n + 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD);
			MPI_Recv(v2, n + 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &status);
			for (j = 0; j < n; j++){
				A[frst * n + j] = v2[j];
			}	
			B[frst] = v2[n];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	BackGauz(A, B, v2, n, my_rank, p, max_rows);
}
