#include "finitevolume.h"
#include "matrix.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main() {

  int N = 8;
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0;
  double courant_fac = 0.4;
  double t = 0.0;
  double tEnd = 2.0;
  double tOut = 0.02;

  double dt = 0.0;
  bool plotThisTurn = false;

  double dx = boxsize / N;
  double vol = pow(dx, 2);

  // Matrixes that live through the entire program -----
  double *rho = (double *)malloc(N * N * sizeof(double));
  double *vx = (double *)malloc(N * N * sizeof(double));
  double *vy = (double *)malloc(N * N * sizeof(double));
  double *P = (double *)malloc(N * N * sizeof(double));

  double *Mass = (double *)malloc(N * N * sizeof(double));
  double *Momx = (double *)malloc(N * N * sizeof(double));
  double *Momy = (double *)malloc(N * N * sizeof(double));
  double *Energy = (double *)malloc(N * N * sizeof(double));

  double *rho_dx = (double *)malloc(N * N * sizeof(double));
  double *rho_dy = (double *)malloc(N * N * sizeof(double));
  double *vx_dx = (double *)malloc(N * N * sizeof(double));
  double *vx_dy = (double *)malloc(N * N * sizeof(double));
  double *vy_dx = (double *)malloc(N * N * sizeof(double));
  double *vy_dy = (double *)malloc(N * N * sizeof(double));
  double *P_dx = (double *)malloc(N * N * sizeof(double));
  double *P_dy = (double *)malloc(N * N * sizeof(double));

  double *rho_prime = (double *)malloc(N * N * sizeof(double));
  double *vx_prime = (double *)malloc(N * N * sizeof(double));
  double *vy_prime = (double *)malloc(N * N * sizeof(double));
  double *P_prime = (double *)malloc(N * N * sizeof(double));
  double *rho_XL = (double *)malloc(N * N * sizeof(double));
  double *rho_XR = (double *)malloc(N * N * sizeof(double));
  double *rho_YL = (double *)malloc(N * N * sizeof(double));
  double *rho_YR = (double *)malloc(N * N * sizeof(double));
  double *vx_XL = (double *)malloc(N * N * sizeof(double));
  double *vx_XR = (double *)malloc(N * N * sizeof(double));
  double *vx_YL = (double *)malloc(N * N * sizeof(double));
  double *vx_YR = (double *)malloc(N * N * sizeof(double));
  double *vy_XL = (double *)malloc(N * N * sizeof(double));
  double *vy_XR = (double *)malloc(N * N * sizeof(double));
  double *vy_YL = (double *)malloc(N * N * sizeof(double));
  double *vy_YR = (double *)malloc(N * N * sizeof(double));
  double *P_XL = (double *)malloc(N * N * sizeof(double));
  double *P_XR = (double *)malloc(N * N * sizeof(double));
  double *P_YL = (double *)malloc(N * N * sizeof(double));
  double *P_YR = (double *)malloc(N * N * sizeof(double));

  double *flux_Mass_X = (double *)malloc(N * N * sizeof(double));
  double *flux_Momx_X = (double *)malloc(N * N * sizeof(double));
  double *flux_Momy_X = (double *)malloc(N * N * sizeof(double));
  double *flux_Energy_X = (double *)malloc(N * N * sizeof(double));

  double *flux_Mass_Y = (double *)malloc(N * N * sizeof(double));
  double *flux_Momx_Y = (double *)malloc(N * N * sizeof(double));
  double *flux_Momy_Y = (double *)malloc(N * N * sizeof(double));
  double *flux_Energy_Y = (double *)malloc(N * N * sizeof(double));
  // ---------------------------------------------------

  // Set initial values for experiment
  double *xlin = (double *)malloc(N * sizeof(double));
  linSpace(0.5 * dx, boxsize - 0.5 * dx, N, xlin);

  double *X = (double *)malloc(N * N * sizeof(double));
  double *Y = (double *)malloc(N * N * sizeof(double));
  meshGrid(xlin, N, X, Y);

  free(xlin);

  double w0 = 0.1;
  double sigma = 0.05 / sqrt(2.0);

  init_rho(Y, N, rho);
  init_vx(Y, N, vx);
  init_vy(X, Y, w0, sigma, N, vy);
  init_P(N, P);

  free(X);
  free(Y);

  getConserved(rho, vx, vy, P, N, gamma, vol, Mass, Momx, Momy, Energy);

  int outputCount = 1;

  while (t < tEnd) {

    getPrimitive(Mass, Momx, Momy, Energy, N, gamma, vol, rho, vx, vy, P);

    dt = getDt(rho, vx, vy, P, N, courant_fac, gamma, dx);

    plotThisTurn = false;

    if (t + dt > outputCount * tOut) {
      dt = outputCount * tOut - t;
      plotThisTurn = true;
    }

    getGradient(rho, N, dx, rho_dx, rho_dy);
    getGradient(vx, N, dx, vx_dx, vx_dy);
    getGradient(vy, N, dx, vy_dx, vy_dy);
    getGradient(P, N, dx, P_dx, P_dy);

    getRhoPrime(rho, rho_dx, rho_dy, vx, vx_dx, vy, vy_dy, dt, N, rho_prime);
    getVxPrime(rho, vx, vx_dx, vx_dy, vy, P_dx, dt, N, vx_prime);
    getVyPrime(rho, vx, vy, vy_dx, vy_dy, P_dy, dt, N, vy_prime);
    getPPrime(vx, vx_dx, vy, vy_dy, P, P_dx, P_dy, dt, gamma, N, P_prime);

    extrapolate(rho_prime, rho_dx, rho_dy, N, dx, rho_XL, rho_XR, rho_YL,
                rho_YR);
    extrapolate(vx_prime, vx_dx, vx_dy, N, dx, vx_XL, vx_XR, vx_YL, vx_YR);
    extrapolate(vy_prime, vy_dx, vy_dy, N, dx, vy_XL, vy_XR, vy_YL, vy_YR);
    extrapolate(P_prime, P_dx, P_dy, N, dx, P_XL, P_XR, P_YL, P_YR);

    getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, N, gamma,
            flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X);
    getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, N, gamma,
            flux_Mass_Y, flux_Momy_Y, flux_Momx_X, flux_Energy_Y);

    applyFluxes(Mass, flux_Mass_X, flux_Mass_Y, N, dx, dt, Mass);
    applyFluxes(Momx, flux_Momx_X, flux_Momx_Y, N, dx, dt, Momx);
    applyFluxes(Momy, flux_Momy_X, flux_Momy_Y, N, dx, dt, Momy);
    applyFluxes(Energy, flux_Energy_X, flux_Energy_Y, N, dx, dt, Energy);

    t += dt;

    if (plotThisTurn || t >= tEnd) {
      outputCount += 1;
      // save or display rho
    }
  }

  printf("Done\n");
  return 0;
}
