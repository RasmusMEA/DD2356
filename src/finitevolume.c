#include "finitevolume.h"
#include "matrix.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void linSpace(double from, double to, int number, double *output) {
  double delta = (to - from) / (number - 1.0);

  int count = 0;
  double curr = from;

  while (count < number) {
    output[count++] = curr;
    curr += delta;
  }
}

void meshGrid(double *arr, int N, double *output1, double *output2) {

  for (int x = 0; x < N; ++x) {
    for (int y = 0; y < N; ++y) {
      output1[y * N + x] = arr[y];
    }
  }

  for (int i = 0; i < N * N; ++i) {
    output2[i] = arr[i % N];
  }
}

void init_rho(double *Y, int N, double *rho_o) {
  copy(rho_o, Y, N);
  sub_MD(rho_o, 0.5, N, rho_o);
  abs_M(rho_o, N, rho_o);
  lessthan_MD(rho_o, 0.25, N, rho_o);
  add_MD(rho_o, 1.0, N, rho_o);
}

void init_vx(double *Y, int N, double *vx_o) {
  copy(vx_o, Y, N);
  sub_MD(vx_o, 0.5, N, vx_o);
  abs_M(vx_o, N, vx_o);
  lessthan_MD(vx_o, 0.25, N, vx_o);
  sub_MD(vx_o, 0.5, N, vx_o);
}

void init_vy(double *X, double *Y, double w0, double sigma, int N,
             double *vy_o) {

  // vy = w0 * np.sin(4 * np.pi * X) *
  //      (np.exp(-(Y - 0.25) * *2 / (2 * sigma * *2)) +
  //       np.exp(-(Y - 0.75) * *2 / (2 * sigma * *2)))

  sub_MD(Y, 0.25, N, vy_o);
  pow2_M(vy_o, N, vy_o);
  div_MD(vy_o, -2.0 * pow(sigma, 2), N, vy_o);
  exp_M(vy_o, N, vy_o);

  double *temp = (double *)malloc(N * N * sizeof(double));

  sub_MD(Y, 0.75, N, temp);
  pow2_M(temp, N, temp);
  div_MD(temp, -2.0 * pow(sigma, 2.0), N, temp);
  exp_M(temp, N, temp);
  add_MM(vy_o, temp, N, vy_o);

  mult_MD(X, 4.0 * M_PI, N, temp);
  sin_M(temp, N, temp);
  mult_MD(temp, w0, N, temp);

  mult_MM(vy_o, temp, N, vy_o);

  free(temp);
}

void init_P(int N, double *P_o) {
  for (int i = 0; i < N * N; ++i) {
    P_o[i] = 2.5;
  }
}

void getConserved(double *rho, double *vx, double *vy, double *P, int N,
                  double gamma, double vol, double *Mass_o, double *Momx_o,
                  double *Momy_o, double *Energy_o) {

  mult_MD(rho, vol, N, Mass_o);

  mult_MM(rho, vx, N, Momx_o);
  mult_MD(Momx_o, vol, N, Momx_o);

  mult_MM(rho, vy, N, Momy_o);
  mult_MD(Momy_o, vol, N, Momy_o);

  double *temp = (double *)malloc(N * N * sizeof(double));

  // calc energy
  pow2_M(vx, N, temp);
  pow2_M(vy, N, Energy_o);

  add_MM(temp, Energy_o, N, temp);
  mult_MM(temp, rho, N, temp);
  mult_MD(temp, 0.5, N, temp);

  div_MD(P, gamma - 1.0, N, Energy_o);
  add_MM(Energy_o, temp, N, Energy_o);
  mult_MD(Energy_o, vol, N, Energy_o);

  free(temp);
}

void getPrimitive(double *Mass, double *Momx, double *Momy, double *Energy,
                  int N, double gamma, double vol, double *rho_o, double *vx_o,
                  double *vy_o, double *P_o) {

  // rho
  div_MD(Mass, vol, N, rho_o);

  // vx
  div_MM(Momx, rho_o, N, vx_o);
  div_MD(vx_o, vol, N, vx_o);

  // vy
  div_MM(Momy, rho_o, N, vy_o);
  div_MD(vy_o, vol, N, vy_o);

  // P
  double *temp = (double *)malloc(N * N * sizeof(double));
  pow2_M(vx_o, N, P_o);
  pow2_M(vy_o, N, temp);
  add_MM(P_o, temp, N, P_o);
  mult_MM(P_o, rho_o, N, P_o);
  mult_MD(P_o, 0.5, N, P_o);

  div_MD(Energy, vol, N, temp);

  sub_MM(temp, P_o, N, P_o);

  mult_MD(P_o, gamma - 1.0, N, P_o);

  free(temp);
}

void getGradient(double *f, int N, double dx, double *f_dx_o, double *f_dy_o) {

  double *temp = (double *)malloc(N * N * sizeof(double));

  rollUp(f, N, f_dx_o);
  rollDown(f, N, temp);
  sub_MM(f_dx_o, temp, N, f_dx_o);
  div_MD(f_dx_o, 2.0 * dx, N, f_dx_o);

  rollLeft(f, N, f_dy_o);
  rollRight(f, N, temp);
  sub_MM(f_dy_o, temp, N, f_dy_o);
  div_MD(f_dy_o, 2.0 * dx, N, f_dy_o);

  free(temp);
}

void getRhoPrime(double *rho, double *rho_dx, double *rho_dy, double *vx,
                 double *vx_dx, double *vy, double *vy_dy, double dt, int N,
                 double *rho_prime_o) {

  // rho_prime = rho - 0.5*dt * (vx * rho_dx + rho * vx_dx + vy * rho_dy + rho *
  // vy_dy)

  double *temp = (double *)malloc(N * N * sizeof(double));
  mult_MM(vx, rho_dx, N, rho_prime_o);
  mult_MM(rho, vx_dx, N, temp);
  add_MM(rho_prime_o, temp, N, rho_prime_o);

  mult_MM(vy, rho_dy, N, temp);
  add_MM(rho_prime_o, temp, N, rho_prime_o);

  mult_MM(rho, vy_dy, N, temp);
  add_MM(rho_prime_o, temp, N, rho_prime_o);

  free(temp);

  mult_MD(rho_prime_o, 0.5 * dt, N, rho_prime_o);

  sub_MM(rho, rho_prime_o, N, rho_prime_o);
}

void getVxPrime(double *rho, double *vx, double *vx_dx, double *vx_dy,
                double *vy, double *P_dx, double dt, int N,
                double *vx_prime_o) {

  double *temp = (double *)malloc(N * N * sizeof(double));

  mult_MM(vx, vx_dx, N, vx_prime_o);
  mult_MM(vy, vx_dy, N, temp);
  add_MM(vx_prime_o, temp, N, vx_prime_o);

  div_DM(1.0, rho, N, temp);
  mult_MM(temp, P_dx, N, temp);
  add_MM(vx_prime_o, temp, N, vx_prime_o);

  mult_MD(vx_prime_o, 0.5 * dt, N, vx_prime_o);

  sub_MM(vx, vx_prime_o, N, vx_prime_o);

  free(temp);
}

void getVyPrime(double *rho, double *vx, double *vy, double *vy_dx,
                double *vy_dy, double *P_dy, double dt, int N,
                double *vy_prime_o) {

  double *temp = (double *)malloc(N * N * sizeof(double));

  mult_MM(vx, vy_dx, N, vy_prime_o);
  mult_MM(vy, vy_dy, N, temp);
  add_MM(vy_prime_o, temp, N, vy_prime_o);

  div_DM(1.0, rho, N, temp);
  mult_MM(temp, P_dy, N, temp);
  add_MM(vy_prime_o, temp, N, vy_prime_o);

  mult_MD(vy_prime_o, 0.5 * dt, N, vy_prime_o);

  sub_MM(vy, vy_prime_o, N, vy_prime_o);

  free(temp);
}

void getPPrime(double *vx, double *vx_dx, double *vy, double *vy_dy, double *P,
               double *P_dx, double *P_dy, double dt, double gamma, int N,
               double *P_prime_o) {

  double *temp = (double *)malloc(N * N * sizeof(double));

  mult_MD(P, gamma, N, P_prime_o);
  add_MM(vx_dx, vy_dy, N, temp);
  mult_MM(P_prime_o, temp, N, P_prime_o);

  mult_MM(vx, P_dx, N, temp);
  add_MM(P_prime_o, temp, N, P_prime_o);

  mult_MM(vy, P_dy, N, temp);
  add_MM(P_prime_o, temp, N, P_prime_o);

  mult_MD(P_prime_o, 0.5 * dt, N, P_prime_o);

  sub_MM(P, P_prime_o, N, P_prime_o);

  free(temp);
}

void extrapolate(double *f, double *f_dx, double *f_dy, int N, double dx,
                 double *f_XL_o, double *f_XR_o, double *f_YL_o,
                 double *f_YR_o) {

  // f_XL
  mult_MD(f_dx, dx / 2.0, N, f_XL_o);
  sub_MM(f, f_XL_o, N, f_XL_o);
  rollUp(f_XL_o, N, f_XL_o);

  // f_XR
  mult_MD(f_dx, dx / 2.0, N, f_XR_o);
  add_MM(f, f_XR_o, N, f_XR_o);

  // f_YL
  mult_MD(f_dy, dx / 2.0, N, f_YL_o);
  sub_MM(f, f_YL_o, N, f_YL_o);
  rollLeft(f_YL_o, N, f_YL_o);

  // f_YR
  mult_MD(f_dy, dx / 2.0, N, f_YR_o);
  add_MM(f, f_YR_o, N, f_YR_o);
}

void getFlux(double *rho_L, double *rho_R, double *vx_L, double *vx_R,
             double *vy_L, double *vy_R, double *P_L, double *P_R, int N,
             double gamma, double *flux_Mass_o, double *flux_Momx_o,
             double *flux_Momy_o, double *flux_Energy_o) {

  double *temp_a = (double *)malloc(N * N * sizeof(double));
  double *temp_b = (double *)malloc(N * N * sizeof(double));

  pow2_M(vx_L, N, temp_a);
  
  // TODO: finish this!!


}

void applyFluxes(double *F, double *flux_F_X, double *flux_F_Y, int N,
                 double dx, double dt, double *output) {

  double *temp = (double *)malloc(N * N * sizeof(double));

  mult_MD(flux_F_X, -dt * dx, N, temp);
  add_MM(F, temp, N, F);

  rollDown(flux_F_X, N, temp);
  mult_MD(temp, dt * dx, N, temp);
  add_MM(F, temp, N, F);

  mult_MD(flux_F_Y, -dt * dx, N, temp);
  add_MM(F, temp, N, F);

  rollRight(flux_F_Y, N, temp);
  mult_MD(temp, dt * dx, N, temp);
  add_MM(F, temp, N, F);

  free(temp);
}

double getDt(double *rho, double *vx, double *vy, double *P, int N,
             double courant_fac, double gamma, double dx) {
  double *temp_a = (double *)malloc(N * N * sizeof(double));
  double *temp_b = (double *)malloc(N * N * sizeof(double));

  pow2_M(vx, N, temp_a);
  pow2_M(vy, N, temp_b);
  add_MM(temp_a, temp_b, N, temp_a);
  sqrt_M(temp_a, N, temp_a);

  mult_MD(P, gamma, N, temp_b);
  div_MM(temp_b, rho, N, temp_b);
  sqrt_M(temp_b, N, temp_b);

  add_MM(temp_a, temp_b, N, temp_a);

  div_DM(dx, temp_a, N, temp_a);

  double new_dt = courant_fac * min(temp_a, N);

  free(temp_a);
  free(temp_b);

  return new_dt;
}

// int main() {
//
//   int N = 7;
//
//   double start = 0.3;
//   double end = 0.7;
//
//   double *xlin = (double *)malloc(N * sizeof(double));
//
//   linSpace(0.2, 0.5, N, xlin);
//
//   double *X = (double *)malloc(N * N * sizeof(double));
//   double *Y = (double *)malloc(N * N * sizeof(double));
//
//   meshGrid(xlin, N, X, Y);
//
//   print(X, N);
//   print(Y, N);
//
//   return 0;
// }
