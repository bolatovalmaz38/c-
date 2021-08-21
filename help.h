#ifndef CREATEMAT_H
#define CREATEMAT_H
double f(int i,int j);
void init_data(double* mat,double* vec,int n,int my_rank,int p,int max_rows, int flag,char *argv[]);
void print_data(double* mat,double* vec,double* v2,int n,int my_rank,int p,int mi,int mj);
#endif
