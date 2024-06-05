#include "finitevolume.h"
#include "matrix.h"
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// debug function for printing a matrix
void printM(double *m, int w, int h) {
  for(int y = 0; y < h; ++y) {
    for(int x = 0; x < w; ++x) {
      printf("%lf ", m[y * w + x]);
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char *argv[]) {
  int rank, size, i, provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int threads = 4;
  omp_set_num_threads(threads);

  int n = 8;

  int N = n;
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

  // TODO: checks to see even divides
  int divWidth = N;
  int divHeight = N / size;

  // allocate local matricies
  double *Mass = (double *)malloc(divWidth * divHeight * sizeof(double));
  double *Momx = (double *)malloc(divWidth * divHeight * sizeof(double));
  double *Momy = (double *)malloc(divWidth * divHeight * sizeof(double));
  double *Energy = (double *)malloc(divWidth * divHeight * sizeof(double));

  // Share initial setup values
  if (rank == 0) {
    // Matricies to hold inital values for experiment
    double *rho_glob = (double *)malloc(N * N * sizeof(double));
    double *vx_glob = (double *)malloc(N * N * sizeof(double));
    double *vy_glob = (double *)malloc(N * N * sizeof(double));
    double *P_glob = (double *)malloc(N * N * sizeof(double));
    double *Mass_glob = (double *)malloc(N * N * sizeof(double));
    double *Momx_glob = (double *)malloc(N * N * sizeof(double));
    double *Momy_glob = (double *)malloc(N * N * sizeof(double));
    double *Energy_glob = (double *)malloc(N * N * sizeof(double));
    double *xlin = (double *)malloc(N * sizeof(double));

    linSpace(0.5 * dx, boxsize - 0.5 * dx, N, xlin);

    double *X = (double *)malloc(N * N * sizeof(double));
    double *Y = (double *)malloc(N * N * sizeof(double));
    meshGrid(xlin, N, X, Y);

    free(xlin);

    double w0 = 0.1;
    double sigma = 0.05 / sqrt(2.0);

    init_rho(Y, N, rho_glob);
    init_vx(Y, N, vx_glob);
    init_vy(X, Y, w0, sigma, N, vy_glob);
    init_P(N, P_glob);

    free(X);
    free(Y);

    getConserved(rho_glob, vx_glob, vy_glob, P_glob, N, gamma, vol, Mass_glob, Momx_glob, Momy_glob, Energy_glob);

    // debug stuff, remove 
    for(int i = 0; i < N * N; ++i) {
      Mass_glob[i] = i + 1;
      Momx_glob[i] = i + 1;
      Momy_glob[i] = i + 1;
      Energy_glob[i] = i + 1;
    }

    printM(Mass_glob, N, N);

    free(rho_glob);
    free(vx_glob);
    free(vy_glob);
    free(P_glob);

    MPI_Request init_requests[4 * (size - 1)];

    //Copy inital values to local array
    memcpy(Mass, Mass_glob, divHeight * divWidth * sizeof(double));
    memcpy(Momx, Momx_glob, divHeight * divWidth * sizeof(double));
    memcpy(Momy, Momy_glob, divHeight * divWidth * sizeof(double));
    memcpy(Energy, Energy_glob, divHeight * divWidth * sizeof(double));

    // send values to all ranks
    for(int i = 1; i < size; ++i) {
      MPI_Isend(&Mass_glob[i * divWidth * divHeight], divHeight * divWidth, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &init_requests[(i - 1) * 4]);
      MPI_Isend(&Momx_glob[i * divWidth * divHeight], divHeight * divWidth, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &init_requests[(i - 1) * 4 + 1]);
      MPI_Isend(&Momy_glob[i * divWidth * divHeight], divHeight * divWidth, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &init_requests[(i - 1) * 4 + 2]);
      MPI_Isend(&Energy_glob[i * divWidth * divHeight], divHeight * divWidth, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &init_requests[(i - 1) * 4 + 3]);
    }


    MPI_Waitall(4 *(size - 1), init_requests, MPI_STATUSES_IGNORE);

    free(Mass_glob);
    free(Momx_glob);
    free(Momy_glob);
    free(Energy_glob);

  } else {
    MPI_Request reqs[4];
    MPI_Irecv(Mass, divHeight * divWidth, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(Momx, divHeight * divWidth, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &reqs[1]);
    MPI_Irecv(Momy, divHeight * divWidth, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(Energy, divHeight * divWidth, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &reqs[3]);

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
  }

  sleep(rank);

  if(rank == 0) {
    printM(Mass, divWidth, divHeight);
  }

  if(rank == 1) {
    printM(Mass, divWidth, divHeight);
  }

  double * rho = (double*) malloc(sizeof(double) * divHeight * divWidth);
  double * vx = (double*) malloc(sizeof(double) * divHeight * divWidth);
  double * vy = (double*) malloc(sizeof(double) * divHeight * divWidth);
  double * P = (double*) malloc(sizeof(double) * divHeight * divWidth);

  while(t < tEnd) {

    getPrimitive(Mass, Momx, Momy, Energy, N, gamma, vol, rho, vx, vy, P);
    dt = getDt(rho, vx, vy, P, N, courant_fac, gamma, dx);
    
    printf("got here\n");

    // sync dt values
    MPI_Send(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    if(rank == 0) {
      double dts[size];
      MPI_Request reqs[size];
      for(int i = 0; i < rank; ++i) {
        MPI_Irecv(&dts[i], size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[i]);
      }
      MPI_Waitall(size, reqs, MPI_STATUSES_IGNORE);
      for(int i = 0; i < size; ++i) {
          dt = fmin(dt, dts[i]);
      }
    }
  }
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  printf("dt: %lf\n", dt);

  // Matrixes that live through the entire program -----
  // double *rho_dx = (double *)malloc(N * N * sizeof(double));
  // double *rho_dy = (double *)malloc(N * N * sizeof(double));
  // double *vx_dx = (double *)malloc(N * N * sizeof(double));
  // double *vx_dy = (double *)malloc(N * N * sizeof(double));
  // double *vy_dx = (double *)malloc(N * N * sizeof(double));
  // double *vy_dy = (double *)malloc(N * N * sizeof(double));
  // double *P_dx = (double *)malloc(N * N * sizeof(double));
  // double *P_dy = (double *)malloc(N * N * sizeof(double));
  // double *rho_prime = (double *)malloc(N * N * sizeof(double));
  // double *vx_prime = (double *)malloc(N * N * sizeof(double));
  // double *vy_prime = (double *)malloc(N * N * sizeof(double));
  // double *P_prime = (double *)malloc(N * N * sizeof(double));
  // double *rho_XL = (double *)malloc(N * N * sizeof(double));
  // double *rho_XR = (double *)malloc(N * N * sizeof(double));
  // double *rho_YL = (double *)malloc(N * N * sizeof(double));
  // double *rho_YR = (double *)malloc(N * N * sizeof(double));
  // double *vx_XL = (double *)malloc(N * N * sizeof(double));
  // double *vx_XR = (double *)malloc(N * N * sizeof(double));
  // double *vx_YL = (double *)malloc(N * N * sizeof(double));
  // double *vx_YR = (double *)malloc(N * N * sizeof(double));
  // double *vy_XL = (double *)malloc(N * N * sizeof(double));
  // double *vy_XR = (double *)malloc(N * N * sizeof(double));
  // double *vy_YL = (double *)malloc(N * N * sizeof(double));
  // double *vy_YR = (double *)malloc(N * N * sizeof(double));
  // double *P_XL = (double *)malloc(N * N * sizeof(double));
  // double *P_XR = (double *)malloc(N * N * sizeof(double));
  // double *P_YL = (double *)malloc(N * N * sizeof(double));
  // double *P_YR = (double *)malloc(N * N * sizeof(double))-;
  // double *flux_Mass_X = (double *)malloc(N * N * sizeof(double));
  // double *flux_Momx_X = (double *)malloc(N * N * sizeof(double));
  // double *flux_Momy_X = (double *)malloc(N * N * sizeof(double));
  // double *flux_Energy_X = (double *)malloc(N * N * sizeof(double));
  // double *flux_Mass_Y = (double *)malloc(N * N * sizeof(double));
  // double *flux_Momx_Y = (double *)malloc(N * N * sizeof(double));
  // double *flux_Momy_Y = (double *)malloc(N * N * sizeof(double));
  // double *flux_Energy_Y = (double *)malloc(N * N * sizeof(double));
  // ---------------------------------------------------

  // getConserved(rho, vx, vy, P, N, gamma, vol, Mass, Momx, Momy, Energy);
  //
  // int outputCount = 1;
  //
  // while (t < tEnd) {
  //
  //   // printf("%lf / %lf\n", t, tEnd);
  //
  //   getPrimitive(Mass, Momx, Momy, Energy, N, gamma, vol, rho, vx, vy, P);
  //
  //   dt = getDt(rho, vx, vy, P, N, courant_fac, gamma, dx);
  //
  //   plotThisTurn = false;
  //
  //   if (t + dt > outputCount * tOut) {
  //     dt = outputCount * tOut - t;
  //     plotThisTurn = true;
  //   }
  //
  //   getGradient(rho, N, dx, rho_dx, rho_dy);
  //   getGradient(vx, N, dx, vx_dx, vx_dy);
  //   getGradient(vy, N, dx, vy_dx, vy_dy);
  //   getGradient(P, N, dx, P_dx, P_dy);
  //
  //   getRhoPrime(rho, rho_dx, rho_dy, vx, vx_dx, vy, vy_dy, dt, N, rho_prime);
  //   getVxPrime(rho, vx, vx_dx, vx_dy, vy, P_dx, dt, N, vx_prime);
  //   getVyPrime(rho, vx, vy, vy_dx, vy_dy, P_dy, dt, N, vy_prime);
  //   getPPrime(vx, vx_dx, vy, vy_dy, P, P_dx, P_dy, dt, gamma, N, P_prime);
  //
  //   extrapolate(rho_prime, rho_dx, rho_dy, N, dx, rho_XL, rho_XR, rho_YL,
  //               rho_YR);
  //   extrapolate(vx_prime, vx_dx, vx_dy, N, dx, vx_XL, vx_XR, vx_YL, vx_YR);
  //   extrapolate(vy_prime, vy_dx, vy_dy, N, dx, vy_XL, vy_XR, vy_YL, vy_YR);
  //   extrapolate(P_prime, P_dx, P_dy, N, dx, P_XL, P_XR, P_YL, P_YR);
  //
  //   getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, N, gamma,
  //           flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X);
  //
  //   getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, N, gamma,
  //           flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y);
  //
  //   applyFluxes(Mass, flux_Mass_X, flux_Mass_Y, N, dx, dt, Mass);
  //   applyFluxes(Momx, flux_Momx_X, flux_Momx_Y, N, dx, dt, Momx);
  //   applyFluxes(Momy, flux_Momy_X, flux_Momy_Y, N, dx, dt, Momy);
  //   applyFluxes(Energy, flux_Energy_X, flux_Energy_Y, N, dx, dt, Energy);
  //
  //   t += dt;
  //
  //   if (plotThisTurn || t >= tEnd) {
  //     outputCount += 1;
  //   }
  // }
  //
  // printf("Done\n");
  return 0;
}
