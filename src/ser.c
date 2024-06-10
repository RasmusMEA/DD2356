#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "finitevol.h"

void simloop(int n) {

  // Simulation parameters
  int N = n; // Resolution of the simulation, i.e x and y dimension of matricies
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0; // Ideal gas gamma
  double courant_fac = 0.4;
  double t = 0.0; // Start time
  double tEnd = 2.0; // End time

  // How often to output frame if real time pltting is desires
  double tOut = 0.02;
  double plotRealTime = true;
  
  // Allocating buffers for all the matricies, these live through the entire duration of the 
  // program and are freed after the simulation loop has finished
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

  // Generate two meshgrids Y and X that are used to
  // generate initial vals for Mass, Momx, Momy and Energy
  meshgrid(xlin, xlin, N, Y, X);

  // Initialize Mass, Momx, Momy and Energy
  for (size_t i = 0; i < N * N; ++i) {
    double val = fabs(Y[i] - 0.5) < 0.25 ? 1.0 : 0.0;

    double rho = val + 1;
    double vx = val - 0.5;
    double vy =
        w0 * sin(4 * M_PI * X[i]) * (exp(-1 * (pow(Y[i] - 0.25, 2) / (2 * pow(sigma, 2)))) + exp(-1 * (pow(Y[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
    double P = 2.5;

    getConserved(rho, vx, vy, P, gamma, vol, &Mass[i], &Momx[i], &Momy[i], &Energy[i]);
  }
  
  // Freeing these since they are only used during
  // initialization
  free(xlin);
  free(X);
  free(Y);

  size_t outputCount = 1;

  // Main simulation loop
  while (t < tEnd) {
    double dt = 999999999999999999.0;

    for (size_t i = 0; i < N * N; ++i) {
      getPrimitive(Mass[i], Momx[i], Momy[i], Energy[i], gamma, vol, &rho[i], &vx[i], &vy[i], &P[i]);

      double dt_i = courant_fac * (dx / (sqrt(gamma * P[i] / rho[i]) + sqrt(pow(vx[i], 2) + pow(vy[i], 2))));

      // Find minimum dt across all values in the matricies
      dt = dt < dt_i ? dt : dt_i;
    }

    // If we should plot this turn, set dt
    // so t is updated to be a multiple of
    // tOut
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

      // print out the final result only
      // change this part of the code if you wish
      // to save every frame of the simulation
      if (t >= tEnd) {
        FILE* stream = fopen("serialoutput.bin", "wb");
        writeToFile(stream, rho, N, N);
        fclose(stream);
      }
    }

    // Calculate gradients and primes for all cells
    for (size_t i = 0; i < N * N; ++i) {
      // Indexing to be able to access
      // elements next to the current element in the matrix
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t down = (y + 1) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;
      size_t right = y * N + (x + 1) % N;

      // Calculate gradients
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

    for (size_t i = 0; i < N * N; ++i) {
      getFlux(rho_XL[i], rho_XR[i], vx_XL[i], vx_XR[i], vy_XL[i], vy_XR[i], P_XL[i], P_XR[i], gamma, &flux_Mass_X[i], &flux_Momx_X[i],
              &flux_Momy_X[i], &flux_Energy_X[i]);

      getFlux(rho_YL[i], rho_YR[i], vy_YL[i], vy_YR[i], vx_YL[i], vx_YR[i], P_YL[i], P_YR[i], gamma, &flux_Mass_Y[i], &flux_Momy_Y[i],
              &flux_Momx_Y[i], &flux_Energy_Y[i]);
    }

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
  
  // Simulation is finished!
  // free all the memory
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

int main(int argc, char* argv[]) {
  int n;
  if (argc < 2) {
    printf("Pass in resulotion as args\n");
    return 0;
  } else {
    n = atoi(argv[1]);
  }

  double t1 = omp_get_wtime();
  simloop(n);
  double t2 = omp_get_wtime();

  printf("time elapsed: %lf seconds\n", t2 - t1);
  return 0;
}
