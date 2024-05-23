#include "matrix.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void _MatrixOnMatrix(double *a, double *b, int N,
                     double (*func)(double, double), double *out) {
  for (size_t i = 0; i < N * N; ++i) {
    out[i] = func(a[i], b[i]);
  }
}

void _MatrixOnDouble(double *a, double b, int N, double (*func)(double, double),
                     double *out) {
  for (size_t i = 0; i < N * N; ++i) {
    out[i] = func(a[i], b);
  }
}

void _MatrixForEach(double *a, int N, double (*func)(double), double *out) {
  for (size_t i = 0; i < N * N; ++i) {
    out[i] = func(a[i]);
  }
}

double add(double a, double b) { return a + b; }

double sub(double a, double b) { return a - b; }

double mult(double a, double b) { return a * b; }

double divi(double a, double b) { return a / b; }

double lessthan(double a, double b) { return (double)(a < b); }

double pow2(double a) { return a * a; }

double expM(double a) { return pow(M_E, a); }

void abs_M(double *a, int N, double *out) { _MatrixForEach(a, N, fabs, out); }

void sin_M(double *a, int N, double *out) { _MatrixForEach(a, N, sin, out); }

void pow2_M(double *a, int N, double *out) { _MatrixForEach(a, N, pow2, out); }

void exp_M(double *a, int N, double *out) { _MatrixForEach(a, N, expM, out); }

void sqrt_M(double *a, int N, double *out) { _MatrixForEach(a, N, sqrt, out); }

void lessthan_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, lessthan, out);
}

void mult_MM(double *a, double *b, int N, double *out) {
  _MatrixOnMatrix(a, b, N, mult, out);
}

void div_MM(double *a, double *b, int N, double *out) {
  _MatrixOnMatrix(a, b, N, divi, out);
}

void add_MM(double *a, double *b, int N, double *out) {
  _MatrixOnMatrix(a, b, N, add, out);
}

void sub_MM(double *a, double *b, int N, double *out) {
  _MatrixOnMatrix(a, b, N, sub, out);
}

void mult_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, mult, out);
}

void div_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, divi, out);
}

void div_DM(double a, double * b, int N, double * out) {
  for(int i = 0; i < N * N; ++i) {
    out[i] = a / b[i];
  }
}

void add_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, add, out);
}

void sub_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, sub, out);
}

void rollUp(double *a, int N, double *out) {
  double *temp = (double *)malloc(N * sizeof(double));

  for (int i = 0; i < N; ++i) {
    temp[i] = a[i];
  }

  for (int y = 0; y < N - 1; ++y) {
    for (int x = 0; x < N; ++x) {
      out[y * N + x] = a[(y + 1) * N + x];
    }
  }

  for (int i = 0; i < N; ++i) {
    out[(N - 1) * N + i] = temp[i];
  }

  free(temp);
}

void rollDown(double *a, int N, double *out) {
  double *temp = (double *)malloc(N * sizeof(double));

  for (int i = 0; i < N; ++i) {
    temp[i] = a[N * (N - 1) + i];
  }

  for (int y = N - 1; y >= 1; --y) {
    for (int x = 0; x < N; ++x) {
      out[y * N + x] = a[(y - 1) * N + x];
    }
  }

  for (int i = 0; i < N; ++i) {
    out[i] = temp[i];
  }

  free(temp);
}

void rollLeft(double *a, int N, double *out) {
  double *temp = (double *)malloc(N * sizeof(double));
  for (int i = 0; i < N; ++i) {
    temp[i] = a[i * N];
  }

  for (int y = 0; y < N; ++y) {
    for (int x = 0; x < N - 1; ++x) {
      out[y * N + x] = a[y * N + x + 1];
    }
  }

  for (int i = 0; i < N; ++i) {
    out[N * i + N - 1] = temp[i];
  }

  free(temp);
}

void rollRight(double *a, int N, double *out) {

  double *temp = (double *)malloc(N * sizeof(double));

  for (int i = 0; i < N; ++i) {
    temp[i] = a[i * N + N - 1];
  }

  for (int y = 0; y < N; ++y) {
    for (int x = N - 1; x >= 1; --x) {
      out[y * N + x] = a[y * N + x - 1];
    }
  }

  for (int i = 0; i < N; ++i) {
    out[N * i] = temp[i];
  }

  free(temp);
}

void print(double *a, int N) {
  for (int y = 0; y < N; ++y) {
    for (int x = 0; x < N; ++x) {
      printf("%lf ", a[y * N + x]);
    }
    printf("\n");
  }
  printf("\n");
}

void copy(double *to, double *from, int N) {
  memcpy(to, from, N * N * sizeof(double));
}

double min(double * a, int N) {
  double min = DBL_MAX;

  for(int i = 0; i < N * N; ++i) {
    if(a[i] < min) {
      min = a[i];
    }
  }

  return min;
}


// int main() {
//
//   int N = 3;
//
//   double *a = (double *)malloc(N * N * sizeof(double));
//
//   double *b = (double *)malloc(N * N * sizeof(double));
//
//   double *c = (double *)malloc(N * N * sizeof(double));
//
//   for (int i = 0; i < N * N; ++i) {
//     a[i] = i + 1.0;
//   }
//
//   printf("a:\n");
//   print(a, N);
//
//   rollUp(a, N, a);
//
//   printf("AFTER:\n");
//
//   printf("a:\n");
//
//   print(a, N);
//
//   usleep(100000);
//
//   free(a);
//   free(b);
//   free(c);
//
//   return 0;
// }
