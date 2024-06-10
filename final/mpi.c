#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "finitevol.h"

// argv[1] = resolution for the simulation
// argv[2] = number of threads for mpi process to use to use
int main(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Pass in Matrix size N and max threads for omp");
  }

  // initalize mpi, only default communicators will be used
  int rank, size, i, provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // set number of processes to be used internally by each node for
  // omp parallelization
  omp_set_num_threads(atoi(argv[2]));

  // variables to measure execution time of the simulation
  double t1, t2;

  if (rank == 0) {
    t1 = omp_get_wtime();
  }

  // Simulation parameters
  int N = atoi(argv[1]);
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0;
  double courant_fac = 0.4;
  double t = 0.0;
  double tEnd = 2.0;
  double tOut = 0.02;
  bool useSlopeLimiting = false;
  double plotRealTime = true;

  int divWidth = N;
  int divHeight = N / size;

  double *rho_glob;
  if (rank == 0) {
    rho_glob = (double *)malloc(sizeof(double) * N * N);
  }

  double *Mass = (double *)malloc(sizeof(double) * divWidth * divHeight);
  double *Momx = (double *)malloc(sizeof(double) * divWidth * divHeight);
  double *Momy = (double *)malloc(sizeof(double) * divWidth * divHeight);
  double *Energy = (double *)malloc(sizeof(double) * divWidth * divHeight);

  // Three extra rows on top and bottom for ghost vals
  double *rho = (double *)malloc(sizeof(double) * divWidth * (divHeight + 6));
  double *vx = (double *)malloc(sizeof(double) * divWidth * (divHeight + 6));
  double *vy = (double *)malloc(sizeof(double) * divWidth * (divHeight + 6));
  double *P = (double *)malloc(sizeof(double) * divWidth * (divHeight + 6));

  // Two extra rows on top and bottom for ghost vals
  double *rho_dx = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *rho_dy = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vx_dx = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vx_dy = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vy_dx = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vy_dy = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *P_dx = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *P_dy = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));

  double *rho_prime = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vx_prime = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vy_prime = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *P_prime = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));

  double *rho_XL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *rho_XR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *rho_YL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *rho_YR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));

  double *vx_XL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vx_XR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vx_YL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vx_YR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));

  double *vy_XL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vy_XR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vy_YL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *vy_YR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));

  double *P_XL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *P_XR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *P_YL = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));
  double *P_YR = (double *)malloc(sizeof(double) * divWidth * (divHeight + 4));

  // One extra row on top and bottom for ghost vals
  double *flux_Mass_X = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));
  double *flux_Momx_X = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));
  double *flux_Momy_X = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));
  double *flux_Energy_X = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));

  double *flux_Mass_Y = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));
  double *flux_Momx_Y = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));
  double *flux_Momy_Y = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));
  double *flux_Energy_Y = (double *)malloc(sizeof(double) * divWidth * (divHeight + 2));

  // Initial conditions
  double w0 = 0.1;
  double sigma = 0.05 / sqrt(2.0);

  // Grid setup
  double dx = boxsize / (double)N;
  double vol = dx * dx;

  // if rank 0, calculate initial setup and distribute
  // values to all mpi processes
  if (rank == 0) {
    double *xlin = (double *)malloc(N * sizeof(double));
    double *X = (double *)malloc(N * N * sizeof(double));
    double *Y = (double *)malloc(N * N * sizeof(double));

    linspace(0.5 * dx, boxsize - 0.5 * dx, N, xlin);

    // Initialize matrices
    meshgrid(xlin, xlin, N, Y, X);

    double *Mass_glob = (double *)malloc(sizeof(double) * N * N);
    double *Momx_glob = (double *)malloc(sizeof(double) * N * N);
    double *Momy_glob = (double *)malloc(sizeof(double) * N * N);
    double *Energy_glob = (double *)malloc(sizeof(double) * N * N);

    // Initialize cells
    for (size_t i = 0; i < N * N; ++i) {
      double val = fabs(Y[i] - 0.5) < 0.25 ? 1.0 : 0.0;

      double rho = val + 1;
      double vx = val - 0.5;
      double vy =
          w0 * sin(4 * M_PI * X[i]) * (exp(-1 * (pow(Y[i] - 0.25, 2) / (2 * pow(sigma, 2)))) + exp(-1 * (pow(Y[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
      double P = 2.5;

      getConserved(rho, vx, vy, P, gamma, vol, &Mass_glob[i], &Momx_glob[i], &Momy_glob[i], &Energy_glob[i]);
    }
    free(xlin);
    free(X);
    free(Y);

    // copy rank 0's initial values directly instead of passing
    // them with mpi
    memcpy(Mass, Mass_glob, sizeof(double) * divWidth * divHeight);
    memcpy(Momx, Momx_glob, sizeof(double) * divWidth * divHeight);
    memcpy(Momy, Momy_glob, sizeof(double) * divWidth * divHeight);
    memcpy(Energy, Energy_glob, sizeof(double) * divWidth * divHeight);

    // Send out initial values to all other ranks
    MPI_Request reqs[4 * (size - 1)];
    for (int i = 1; i < size; ++i) {
      MPI_Isend(&Mass_glob[i * divWidth * divHeight], divWidth * divHeight, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[(i - 1) * 4]);
      MPI_Isend(&Momx_glob[i * divWidth * divHeight], divWidth * divHeight, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &reqs[(i - 1) * 4 + 1]);
      MPI_Isend(&Momy_glob[i * divWidth * divHeight], divWidth * divHeight, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &reqs[(i - 1) * 4 + 2]);
      MPI_Isend(&Energy_glob[i * divWidth * divHeight], divWidth * divHeight, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &reqs[(i - 1) * 4 + 3]);
    }

    MPI_Waitall(4 * (size - 1), reqs, MPI_STATUSES_IGNORE);

    free(Mass_glob);
    free(Momx_glob);
    free(Momy_glob);
    free(Energy_glob);

  } else {
    // if not rank 0, wait to recieve all inital values
    MPI_Request reqs[4];

    MPI_Irecv(Mass, divWidth * divHeight, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(Momx, divWidth * divHeight, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &reqs[1]);
    MPI_Irecv(Momy, divWidth * divHeight, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(Energy, divWidth * divHeight, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &reqs[3]);

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
  }

  size_t outputCount = 1;

  // Main simulation loop
  while (t < tEnd) {
    double dt = 999999999999999999.0;
    double *dt_i = (double *)malloc(sizeof(double) * omp_get_max_threads());

    #pragma omp parallel for
    for (int i = 0; i < omp_get_max_threads(); ++i) {
      dt_i[i] = dt;
    }

    #pragma omp parallel for
    for (size_t i = 0; i < divWidth * divHeight; ++i) {
      // Get primitive values
      getPrimitive(Mass[i], Momx[i], Momy[i], Energy[i], gamma, vol, &rho[3 * divWidth + i], &vx[3 * divWidth + i], &vy[3 * divWidth + i],
                   &P[divWidth * 3 + i]);

      // calculate dt to determine timestep later
      double dt =
          courant_fac *
          (dx / (sqrt(gamma * P[divWidth * 3 + i] / rho[divWidth * 3 + i]) + sqrt(pow(vx[divWidth * 3 + i], 2) + pow(vy[divWidth * 3 + i], 2))));

      if (dt < dt_i[omp_get_thread_num()]) {
        dt_i[omp_get_thread_num()] = dt;
      }
    }

    // Calculate local dt minimum for this mpi process
    dt = dt_i[0];
    for (int i = 1; i < omp_get_max_threads(); ++i) {
      if (dt_i[i] < dt) {
        dt = dt_i[i];
      }
    }
    free(dt_i);

    // Sync minumum dt across all mpi nodes

    // if rank 0, recieve all mpi processes local min dt's and respond with the global min
    if (rank == 0) {
      double dt_i;
      for (int i = 1; i < size; ++i) {
        MPI_Recv(&dt_i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (dt_i < dt) {
          dt = dt_i;
        }
      }
    } else {
      // if not rank 0, send your local min dt to rank 0
      MPI_Send(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    // Send and recieve global min dt
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
      // If we need to save vals this turn, rank 0 will accumulate all
      // local rho values
      if (t >= tEnd) {
        if (rank == 0) {
          // copy rank 0 values directly instead of sending
          memcpy(rho_glob, &rho[divWidth * 3], divWidth * divHeight * sizeof(double));

          // recieve rho values from all other ranks and place them
          // in the correct spot of the global matrix
          MPI_Request reqs[size - 1];
          for (int i = 1; i < size; ++i) {
            MPI_Irecv(&rho_glob[i * divWidth * divHeight], divWidth * divHeight, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &reqs[i - 1]);
          }

          MPI_Waitall(size - 1, reqs, MPI_STATUSES_IGNORE);

          // Write the result in binary to a file
          FILE *stream = fopen("mpioutput.bin", "wb");
          writeToFile(stream, rho_glob, N, N);
          fclose(stream);

        } else {
          // if not rank 0, send local rho vals.
          // offset by divWidth * 3 since the first three rows in rho are
          // ghost values
          MPI_Send(&rho[divWidth * 3], divWidth * divHeight, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
      }

      outputCount += 1;
    }

    // Send and recieve ghost values for rho, vx, vy and P to/from MPI processes
    // above and below our rank.

    // Send three top rows of rho, vx, vy and P to rank above
    MPI_Request reqs[16];
    MPI_Isend(&rho[divWidth * 3], divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &reqs[0]);
    MPI_Isend(&vx[divWidth * 3], divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(&vy[divWidth * 3], divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 2, MPI_COMM_WORLD, &reqs[2]);
    MPI_Isend(&P[divWidth * 3], divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 3, MPI_COMM_WORLD, &reqs[3]);

    // Send three bottom rows of rho, vx, vy and P to rank below
    MPI_Isend(&rho[divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 4, MPI_COMM_WORLD, &reqs[4]);
    MPI_Isend(&vx[divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 5, MPI_COMM_WORLD, &reqs[5]);
    MPI_Isend(&vy[divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 6, MPI_COMM_WORLD, &reqs[6]);
    MPI_Isend(&P[divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 7, MPI_COMM_WORLD, &reqs[7]);

    // Recieving ghost vals from above
    MPI_Irecv(rho, divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 4, MPI_COMM_WORLD, &reqs[8]);
    MPI_Irecv(vx, divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 5, MPI_COMM_WORLD, &reqs[9]);
    MPI_Irecv(vy, divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 6, MPI_COMM_WORLD, &reqs[10]);
    MPI_Irecv(P, divWidth * 3, MPI_DOUBLE, (rank - 1 + size) % size, 7, MPI_COMM_WORLD, &reqs[11]);

    // Recieving ghost vals from below
    MPI_Irecv(&rho[divWidth * 3 + divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &reqs[12]);
    MPI_Irecv(&vx[divWidth * 3 + divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 1, MPI_COMM_WORLD, &reqs[13]);
    MPI_Irecv(&vy[divWidth * 3 + divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 2, MPI_COMM_WORLD, &reqs[14]);
    MPI_Irecv(&P[divWidth * 3 + divWidth * divHeight], divWidth * 3, MPI_DOUBLE, (rank + 1) % size, 3, MPI_COMM_WORLD, &reqs[15]);

    MPI_Waitall(16, reqs, MPI_STATUSES_IGNORE);

    // Skipping first and last row since we don't need to calculate these values,
    // we only need them to be in the matrix because the rest of the matrix
    // depends on them
    #pragma omp parallel for
    for (int y = 1; y < divHeight + 5; ++y) {
      for (int x = 0; x < divWidth; ++x) {
        int up = ((y - 1 + divHeight + 6) % (divHeight + 6)) * divWidth + x;
        int down = ((y + 1) % (divHeight + 6)) * divWidth + x;
        int left = y * divWidth + ((x - 1 + divWidth) % divWidth);
        int right = y * divWidth + ((x + 1) % divWidth);
        int i = y * divWidth + x;

        getGradient(rho[i], rho[up], rho[down], rho[left], rho[right], dx, &rho_dx[i - divWidth], &rho_dy[i - divWidth]);
        getGradient(vx[i], vx[up], vx[down], vx[left], vx[right], dx, &vx_dx[i - divWidth], &vx_dy[i - divWidth]);
        getGradient(vy[i], vy[up], vy[down], vy[left], vy[right], dx, &vy_dx[i - divWidth], &vy_dy[i - divWidth]);
        getGradient(P[i], P[up], P[down], P[left], P[right], dx, &P_dx[i - divWidth], &P_dy[i - divWidth]);

        rho_prime[i - divWidth] = rho[i] - 0.5 * dt *
                                               ((vx[i] * rho_dx[i - divWidth]) + (rho[i] * vx_dx[i - divWidth]) + (vy[i] * rho_dy[i - divWidth]) +
                                                (rho[i] * vy_dy[i - divWidth]));

        vx_prime[i - divWidth] = vx[i] - (0.5 * dt) * (vx[i] * vx_dx[i - divWidth] + vy[i] * vx_dy[i - divWidth] + (1 / rho[i]) * P_dx[i - divWidth]);
        vy_prime[i - divWidth] = vy[i] - (0.5 * dt) * (vx[i] * vy_dx[i - divWidth] + vy[i] * vy_dy[i - divWidth] + (1 / rho[i]) * P_dy[i - divWidth]);
        P_prime[i - divWidth] = P[i] - (0.5 * dt) * (gamma * P[i] * (vx_dx[i - divWidth] + vy_dy[i - divWidth]) + vx[i] * P_dx[i - divWidth] +
                                                     vy[i] * P_dy[i - divWidth]);
      }
    }

    #pragma omp parallel for
    for (int y = 0; y < divHeight + 4; ++y) {
      for (int x = 0; x < divWidth; ++x) {
        int up = ((y - 1 + divHeight + 4) % (divHeight + 4)) * divWidth + x;
        int left = y * divWidth + ((x - 1 + divWidth) % divWidth);
        int i = y * divWidth + x;

        extrapolateInSpaceToFace(rho_prime[i], rho_dx[i], rho_dy[i], dx, &rho_XL[up], &rho_XR[i], &rho_YL[left], &rho_YR[i]);
        extrapolateInSpaceToFace(vx_prime[i], vx_dx[i], vx_dy[i], dx, &vx_XL[up], &vx_XR[i], &vx_YL[left], &vx_YR[i]);
        extrapolateInSpaceToFace(vy_prime[i], vy_dx[i], vy_dy[i], dx, &vy_XL[up], &vy_XR[i], &vy_YL[left], &vy_YR[i]);
        extrapolateInSpaceToFace(P_prime[i], P_dx[i], P_dy[i], dx, &P_XL[up], &P_XR[i], &P_YL[left], &P_YR[i]);
      }
    }

    #pragma omp parallel for
    for (int y = 1; y < divHeight + 3; ++y) {
      for (int x = 0; x < divWidth; ++x) {
        int i = y * divWidth + x;

        getFlux(rho_XL[i], rho_XR[i], vx_XL[i], vx_XR[i], vy_XL[i], vy_XR[i], P_XL[i], P_XR[i], gamma, &flux_Mass_X[i - divWidth],
                &flux_Momx_X[i - divWidth], &flux_Momy_X[i - divWidth], &flux_Energy_X[i - divWidth]);

        getFlux(rho_YL[i], rho_YR[i], vy_YL[i], vy_YR[i], vx_YL[i], vx_YR[i], P_YL[i], P_YR[i], gamma, &flux_Mass_Y[i - divWidth],
                &flux_Momy_Y[i - divWidth], &flux_Momx_Y[i - divWidth], &flux_Energy_Y[i - divWidth]);
      }
    }
    
    #pragma omp parallel for
    for (int y = 1; y < divHeight + 1; ++y) {
      for (int x = 0; x < divWidth; ++x) {
        int i = y * divWidth + x;
        int up = ((y - 1 + divHeight + 2) % (divHeight + 2)) * divWidth + x;
        int left = y * divWidth + ((x - 1 + divWidth) % (divWidth));

        Mass[i - divWidth] = applyFluxes(Mass[i - divWidth], flux_Mass_X[i], flux_Mass_X[up], flux_Mass_Y[left], flux_Mass_Y[i], dx, dt);
        Momx[i - divWidth] = applyFluxes(Momx[i - divWidth], flux_Momx_X[i], flux_Momx_X[up], flux_Momx_Y[left], flux_Momx_Y[i], dx, dt);
        Momy[i - divWidth] = applyFluxes(Momy[i - divWidth], flux_Momy_X[i], flux_Momy_X[up], flux_Momy_Y[left], flux_Momy_Y[i], dx, dt);
        Energy[i - divWidth] = applyFluxes(Energy[i - divWidth], flux_Energy_X[i], flux_Energy_X[up], flux_Energy_Y[left], flux_Energy_Y[i], dx, dt);
      }
    }
  }

  if (rank == 0) {
    t2 = omp_get_wtime();
    double res = t2 - t1;
    printf("execution time: %lf\n", res);
  }

  if (rank == 0) {
    free(rho_glob);
  }
  free(Mass);
  free(Momx);
  free(Momy);
  free(Energy);
  free(rho);
  free(vx);
  free(vy);
  free(P);
  free(rho_dx);
  free(rho_dy);
  free(vx_dx);
  free(vx_dy);
  free(vy_dx);
  free(vy_dy);
  free(P_dx);
  free(P_dy);
  free(rho_prime);
  free(vx_prime);
  free(vy_prime);
  free(P_prime);
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
  free(flux_Momx_X);
  free(flux_Momy_X);
  free(flux_Energy_X);
  free(flux_Mass_Y);
  free(flux_Momx_Y);
  free(flux_Momy_Y);
  free(flux_Energy_Y);

  MPI_Finalize();
  return 0;
}
