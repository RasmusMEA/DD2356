#include "matrix.h"
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void _MatrixOnMatrix(double *a, double *b, int N,
                     double (*func)(double, double), double *out) {

#pragma omp parallel for
  for (size_t i = 0; i < N; ++i) {
    out[i] = func(a[i], b[i]);
  }
}

void _MatrixOnDouble(double *a, double b, int N, double (*func)(double, double),
                     double *out) {
#pragma omp parallel for
  for (size_t i = 0; i < N; ++i) {
    out[i] = func(a[i], b);
  }
}

void _MatrixForEach(double *a, int N, double (*func)(double), double *out) {
#pragma omp parallel for
  for (size_t i = 0; i <N; ++i) {
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

double max(double a, double b) { return a > b ? a : b; }

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

void max_MM(double *a, double *b, int N, double *out) {
  _MatrixOnMatrix(a, b, N, max, out);
}

void mult_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, mult, out);
}

void div_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, divi, out);
}

void div_DM(double a, double *b, int N, double *out) {
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    out[i] = a / b[i];
  }
}

void add_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, add, out);
}

void sub_MD(double *a, double b, int N, double *out) {
  _MatrixOnDouble(a, b, N, sub, out);
}

void rollUp(double *a, int w, int h, double *out) {
  double *temp = (double *)malloc(w * sizeof(double));

#pragma omp parallel for
  for (int i = 0; i < w; ++i) {
    temp[i] = a[i];
  }

#pragma omp parallel for
  for (int y = 0; y < h - 1; ++y) {
    for (int x = 0; x < w; ++x) {
      out[y * w + x] = a[(y + 1) * w + x];
    }
  }

#pragma omp parallel for
  for (int i = 0; i < w; ++i) {
    out[(h - 1) * w + i] = temp[i];
  }

  free(temp);
}

void rollDown(double *a, int w, int h, double *out) {
  double *temp = (double *)malloc(w * sizeof(double));

#pragma omp parallel for
  for (int i = 0; i < w; ++i) {
    temp[i] = a[h * (w - 1) + i];
  }

#pragma omp parallel for
  for (int y = h - 1; y >= 1; --y) {
    for (int x = 0; x < w; ++x) {
      out[y * w + x] = a[(y - 1) * w + x];
    }
  }

#pragma omp parallel for
  for (int i = 0; i < w; ++i) {
    out[i] = temp[i];
  }

  free(temp);
}

void rollLeft(double *a, int w, int h, double *out) {
  double *temp = (double *)malloc(h * sizeof(double));
#pragma omp parallel for
  for (int i = 0; i < h; ++i) {
    temp[i] = a[i * w];
  }

#pragma omp parallel for
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < h - 1; ++x) {
      out[y * w + x] = a[y * w + x + 1];
    }
  }

#pragma omp parallel for
  for (int i = 0; i < h; ++i) {
    out[h * i + w - 1] = temp[i];
  }

  free(temp);
}

void rollRight(double *a, int w, int h, double *out) {

  double *temp = (double *)malloc(h * sizeof(double));

#pragma omp parallel for
  for (int i = 0; i < h; ++i) {
    temp[i] = a[i * h + w - 1];
  }

#pragma omp parallel for
  for (int y = 0; y < h; ++y) {
    for (int x = w - 1; x >= 1; --x) {
      out[y * w + x] = a[y * w + x - 1];
    }
  }

#pragma omp parallel for
  for (int i = 0; i < h; ++i) {
    out[w * i] = temp[i];
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

double min(double *a, int N) {
  double min = DBL_MAX;

#pragma omp parallel for
  for (int i = 0; i < N * N; ++i) {
    if (a[i] < min) {
      min = a[i];
    }
  }

  return min;
}
