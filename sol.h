#include"mspin.h"
#ifndef MATRMUL_H
#define MATRMUL_H
void Norma(double *a, double *b, int n, double *x, int r, int row);
int spin(double *v1, double *v2, int k, mspin *s);
void rotate(double *v1, double *v2, double *w1, double *w2, int n, int k, mspin *ms);
void Solve_it(double *A, double *B,double *v2, int n, mspin *ms, int my_rank, int p, int max_rows);
#endif
