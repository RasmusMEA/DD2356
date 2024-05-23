#ifndef MATRIX_H
#define MATRIX_H

void abs_M(double *a, int N, double *out);
void sin_M(double *a, int N, double *out);
void pow2_M(double *a, int N, double *out);
void exp_M(double *a, int N, double *out);
void sqrt_M(double *a, int N, double *out);

void lessthan_MD(double *a, double b, int N, double *out);

void mult_MM(double *a, double *b, int N, double *out);
void div_MM(double *a, double *b, int N, double *out);
void add_MM(double *a, double *b, int N, double *out);
void sub_MM(double *a, double *b, int N, double *out);

void mult_MD(double *a, double b, int N, double *out);
void div_MD(double *a, double b, int N, double *out);
void div_DM(double a, double * b, int N, double * out);
void add_MD(double *a, double b, int N, double *out);
void sub_MD(double *a, double b, int N, double *out);

void rollUp(double *a, int N, double *out);
void rollDown(double *a, int N, double *out);
void rollLeft(double *a, int N, double *out);
void rollRight(double *a, int N, double *out);

double min(double * a, int N);

void copy(double *to, double *from, int N);
void print(double *a, int N);

#endif
