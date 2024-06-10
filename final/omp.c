#include "finitevol.h"
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void simloop(int n) {
  // Simulation parameters
  int N = n;
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0;
  double courant_fac = 0.4;
  double t = 0.0;
  double tEnd = 2.0;
  double tOut = 0.02;
  double plotRealTime = true;

  double* Mass = (double*)malloc(sizeof(double) * N * N);
  double* Momx = (double*)malloc(sizeof(double) * N * N);
  double* Momy = (double*)malloc(sizeof(double) * N * N);
  double* Energy = (double*)malloc(sizeof(double) * N * N);
  double* rho = (double*)malloc(sizeof(double) * N * N);
  double* vx = (double*)malloc(sizeof(double) * N * N);
  double* vy = (double*)malloc(sizeof(double) * N * N);
  double* P = (double*)malloc(sizeof(double) * N * N);
  double* rho_prime = (double*)malloc(sizeof(double) * N * N);
  double* vx_prime = (double*)malloc(sizeof(double) * N * N);
  double* vy_prime = (double*)malloc(sizeof(double) * N * N);
  double* P_prime = (double*)malloc(sizeof(double) * N * N);
  double* rho_dx = (double*)malloc(sizeof(double) * N * N);
  double* rho_dy = (double*)malloc(sizeof(double) * N * N);
  double* vx_dx = (double*)malloc(sizeof(double) * N * N);
  double* vx_dy = (double*)malloc(sizeof(double) * N * N);
  double* vy_dx = (double*)malloc(sizeof(double) * N * N);
  double* vy_dy = (double*)malloc(sizeof(double) * N * N);
  double* P_dx = (double*)malloc(sizeof(double) * N * N);
  double* P_dy = (double*)malloc(sizeof(double) * N * N);
  double* rho_XL = (double*)malloc(sizeof(double) * N * N);
  double* rho_XR = (double*)malloc(sizeof(double) * N * N);
  double* rho_YL = (double*)malloc(sizeof(double) * N * N);
  double* rho_YR = (double*)malloc(sizeof(double) * N * N);
  double* vx_XL = (double*)malloc(sizeof(double) * N * N);
  double* vx_XR = (double*)malloc(sizeof(double) * N * N);
  double* vx_YL = (double*)malloc(sizeof(double) * N * N);
  double* vx_YR = (double*)malloc(sizeof(double) * N * N);
  double* vy_XL = (double*)malloc(sizeof(double) * N * N);
  double* vy_XR = (double*)malloc(sizeof(double) * N * N);
  double* vy_YL = (double*)malloc(sizeof(double) * N * N);
  double* vy_YR = (double*)malloc(sizeof(double) * N * N);
  double* P_XL = (double*)malloc(sizeof(double) * N * N);
  double* P_XR = (double*)malloc(sizeof(double) * N * N);
  double* P_YL = (double*)malloc(sizeof(double) * N * N);
  double* P_YR = (double*)malloc(sizeof(double) * N * N);
  double* flux_Mass_X = (double*)malloc(sizeof(double) * N * N);
  double* flux_Mass_Y = (double*)malloc(sizeof(double) * N * N);
  double* flux_Momx_X = (double*)malloc(sizeof(double) * N * N);
  double* flux_Momx_Y = (double*)malloc(sizeof(double) * N * N);
  double* flux_Momy_X = (double*)malloc(sizeof(double) * N * N);
  double* flux_Momy_Y = (double*)malloc(sizeof(double) * N * N);
  double* flux_Energy_X = (double*)malloc(sizeof(double) * N * N);
  double* flux_Energy_Y = (double*)malloc(sizeof(double) * N * N);

  // Initial conditions
  double w0 = 0.1;
  double sigma = 0.05 / sqrt(2.0);

  // Grid setup
  double dx = boxsize / (double)N;
  double vol = dx * dx;

  double* xlin = (double*)malloc(N * sizeof(double));
  double* X = (double*)malloc(N * N * sizeof(double));
  double* Y = (double*)malloc(N * N * sizeof(double));

  linspace(0.5 * dx, boxsize - 0.5 * dx, N, xlin);

  // Initialize matrices
  meshgrid(xlin, xlin, N, Y, X);

  // Initialize cells
  for (size_t i = 0; i < N * N; ++i) {
    double val = fabs(Y[i] - 0.5) < 0.25 ? 1.0 : 0.0;

    double rho = val + 1;
    double vx = val - 0.5;
    double vy =
        w0 * sin(4 * M_PI * X[i]) * (exp(-1 * (pow(Y[i] - 0.25, 2) / (2 * pow(sigma, 2)))) + exp(-1 * (pow(Y[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
    double P = 2.5;

    getConserved(rho, vx, vy, P, gamma, vol, &Mass[i], &Momx[i], &Momy[i], &Energy[i]);
  }

  free(xlin);
  free(X);
  free(Y);

  size_t outputCount = 1;

  // Main simulation loop
  while (t < tEnd) {
    double dt = 999999999999999999.0;

    // buffer to save local minimum for dt to avoid
    // using a critical section in the getPrimitive loop
    double* dt_i = (double*)malloc(sizeof(double) * omp_get_max_threads());

    // initalize all values in dti to be a large double, so we can find local minimum
    #pragma omp parallel for
    for (int i = 0; i < omp_get_max_threads(); ++i) {
      dt_i[i] = dt;
    }

    
    #pragma omp parallel for
    for (size_t i = 0; i < N * N; ++i) {
      // Calculate primitive variables
      getPrimitive(Mass[i], Momx[i], Momy[i], Energy[i], gamma, vol, &rho[i], &vx[i], &vy[i], &P[i]);

      // calculate timestep in order to determine
      // minimum timestep later
      double dt = courant_fac * (dx / (sqrt(gamma * P[i] / rho[i]) + sqrt(pow(vx[i], 2) + pow(vy[i], 2))));

      if (dt < dt_i[omp_get_thread_num()]) {
        dt_i[omp_get_thread_num()] = dt;
      }
    }

    // Find global minimum dt across all processes
    dt = dt_i[0];
    for (int i = 1; i < omp_get_max_threads(); ++i) {
      if (dt_i[i] < dt) {
        dt = dt_i[i];
      }
    }
    free(dt_i);

    // Check if we should plot this turn
    bool plotThisTurn = false;
    if (t + dt > outputCount * tOut) {
      dt = outputCount * tOut - t;
      plotThisTurn = true;
    }

    // Update time
    t += dt;

    // Plot
    if (plotThisTurn || t >= tEnd) {
      outputCount += 1;

      if( t >= tEnd) {
        FILE * stream = fopen("ompoutput.bin", "wb");
        writeToFile(stream, rho, N, N);
        fclose(stream);
      }

    }

    // Calculate next step for all cells
    #pragma omp parallel for
    for (size_t i = 0; i < N * N; ++i) {
      // Indexing
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t down = (y + 1) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;
      size_t right = y * N + (x + 1) % N;

      // Requires nearby cells
      getGradient(rho[i], rho[up], rho[down], rho[left], rho[right], dx, &rho_dx[i], &rho_dy[i]);
      getGradient(vx[i], vx[up], vx[down], vx[left], vx[right], dx, &vx_dx[i], &vx_dy[i]);
      getGradient(vy[i], vy[up], vy[down], vy[left], vy[right], dx, &vy_dx[i], &vy_dy[i]);
      getGradient(P[i], P[up], P[down], P[left], P[right], dx, &P_dx[i], &P_dy[i]);

      // Extrapolate half step in time
      rho_prime[i] = rho[i] - 0.5 * dt * ((vx[i] * rho_dx[i]) + (rho[i] * vx_dx[i]) + (vy[i] * rho_dy[i]) + (rho[i] * vy_dy[i]));

      vx_prime[i] = vx[i] - (0.5 * dt) * (vx[i] * vx_dx[i] + vy[i] * vx_dy[i] + (1 / rho[i]) * P_dx[i]);
      vy_prime[i] = vy[i] - (0.5 * dt) * (vx[i] * vy_dx[i] + vy[i] * vy_dy[i] + (1 / rho[i]) * P_dy[i]);
      P_prime[i] = P[i] - (0.5 * dt) * (gamma * P[i] * (vx_dx[i] + vy_dy[i]) + vx[i] * P_dx[i] + vy[i] * P_dy[i]);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < N * N; ++i) {
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;

      extrapolateInSpaceToFace(rho_prime[i], rho_dx[i], rho_dy[i], dx, &rho_XL[up], &rho_XR[i], &rho_YL[left], &rho_YR[i]);
      extrapolateInSpaceToFace(vx_prime[i], vx_dx[i], vx_dy[i], dx, &vx_XL[up], &vx_XR[i], &vx_YL[left], &vx_YR[i]);
      extrapolateInSpaceToFace(vy_prime[i], vy_dx[i], vy_dy[i], dx, &vy_XL[up], &vy_XR[i], &vy_YL[left], &vy_YR[i]);
      extrapolateInSpaceToFace(P_prime[i], P_dx[i], P_dy[i], dx, &P_XL[up], &P_XR[i], &P_YL[left], &P_YR[i]);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < N * N; ++i) {
      getFlux(rho_XL[i], rho_XR[i], vx_XL[i], vx_XR[i], vy_XL[i], vy_XR[i], P_XL[i], P_XR[i], gamma, &flux_Mass_X[i], &flux_Momx_X[i],
              &flux_Momy_X[i], &flux_Energy_X[i]);

      getFlux(rho_YL[i], rho_YR[i], vy_YL[i], vy_YR[i], vx_YL[i], vx_YR[i], P_YL[i], P_YR[i], gamma, &flux_Mass_Y[i], &flux_Momy_Y[i],
              &flux_Momx_Y[i], &flux_Energy_Y[i]);
    }

    #pragma omp parallel for
    for (size_t i = 0; i < N * N; ++i) {
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;

      Mass[i] = applyFluxes(Mass[i], flux_Mass_X[i], flux_Mass_X[up], flux_Mass_Y[left], flux_Mass_Y[i], dx, dt);
      Momx[i] = applyFluxes(Momx[i], flux_Momx_X[i], flux_Momx_X[up], flux_Momx_Y[left], flux_Momx_Y[i], dx, dt);
      Momy[i] = applyFluxes(Momy[i], flux_Momy_X[i], flux_Momy_X[up], flux_Momy_Y[left], flux_Momy_Y[i], dx, dt);
      Energy[i] = applyFluxes(Energy[i], flux_Energy_X[i], flux_Energy_X[up], flux_Energy_Y[left], flux_Energy_Y[i], dx, dt);
    }
  }

  free(Mass);
  free(Momx);
  free(Momy);
  free(Energy);
  free(rho);
  free(vx);
  free(vy);
  free(P);
  free(rho_prime);
  free(vx_prime);
  free(vy_prime);
  free(P_prime);
  free(rho_dx);
  free(rho_dy);
  free(vx_dx);
  free(vx_dy);
  free(vy_dx);
  free(vy_dy);
  free(P_dx);
  free(P_dy);
  free(rho_XL);
  free(rho_XR);
  free(rho_YL);
  free(rho_YR);
  free(vx_XL);
  free(vx_XR);
  free(vx_YL);
  free(vx_YR);
  free(vy_XL);
  free(vy_XR);
  free(vy_YL);
  free(vy_YR);
  free(P_XL);
  free(P_XR);
  free(P_YL);
  free(P_YR);
  free(flux_Mass_X);
  free(flux_Mass_Y);
  free(flux_Momx_X);
  free(flux_Momx_Y);
  free(flux_Momy_X);
  free(flux_Momy_Y);
  free(flux_Energy_X);
  free(flux_Energy_Y);
}

// argv[1] = resolution for the simulation
// argv[2] = number of threads for omp to use
int main(int argc, char* argv[]) {
  int n;
  if (argc < 3) {
    printf("Pass in resulotion as args\n");
    return 0;
  }

  omp_set_num_threads(atoi(argv[2]));

  double t1 = omp_get_wtime();
  simloop(atoi(argv[1]));
  double t2 = omp_get_wtime();

  printf("time elapsed: %lf seconds\n", t2 - t1);
  return 0;
}
