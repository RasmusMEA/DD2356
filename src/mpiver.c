#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


void getConserved(double rho, double vx, double vy, double P, double gamma,
                  double vol, double *Mass_o, double *Momx_o, double *Momy_o,
                  double *Energy_o) {
  double Mass = rho * vol;
  double Momx = Mass * vx;
  double Momy = Mass * vy;
  double Energy = (P / (gamma - 1) + 0.5 * rho * ((vx * vx) + (vy * vy))) * vol;

  *Mass_o = Mass;
  *Momx_o = Momx;
  *Momy_o = Momy;
  *Energy_o = Energy;
}

void getPrimitive(double Mass, double Momx, double Momy, double Energy,
                  double gamma, double vol, double *rho_o, double *vx_o,
                  double *vy_o, double *P_o) {
  double rho = Mass / vol;
  double vx = Momx / rho / vol;
  double vy = Momy / rho / vol;
  double P = (Energy / vol - 0.5 * rho * ((vx * vx) + (vy * vy))) * (gamma - 1);

  *rho_o = rho;
  *vx_o = vx;
  *vy_o = vy;
  *P_o = P;
}

void getGradient(double f_down, double f_up, double f_right,
                 double f_left, double dx, double *f_dx_o, double *f_dy_o) {
  double f_dx = (f_up - f_down) / (2 * dx);
  double f_dy = (f_left - f_right) / (2 * dx);

  *f_dx_o = f_dx;
  *f_dy_o = f_dy;
}

void extrapolateInSpaceToFace(double f, double f_dx, double f_dy, double dx,
                              double *f_XL_o, double *f_XR_o, double *f_YL_o,
                              double *f_YR_o) {
  double f_XL = f - f_dx * (dx / 2.0);
  double f_XR = f + f_dx * (dx / 2.0);

  double f_YL = f - f_dy * (dx / 2.0);
  double f_YR = f + f_dy * (dx / 2.0);

  *f_XL_o = f_XL;
  *f_XR_o = f_XR;
  *f_YR_o = f_YR;
  *f_YL_o = f_YL;
}

double applyFluxes(double F, double flux_F_x, double flux_F_x_up,
                   double flux_F_y_left, double flux_F_y, double dx,
                   double dt) {
  return F + (dt * dx) * (-flux_F_x + flux_F_x_up - flux_F_y + flux_F_y_left);
}

void getFlux(double rho_L, double rho_R, double vx_L, double vx_R, double vy_L,
             double vy_R, double P_L, double P_R, double gamma,
             double *flux_Mass_o, double *flux_Momx_o, double *flux_Momy_o,
             double *flux_Energy_o) {
  double en_L =
      P_L / (gamma - 1) + 0.5 * rho_L * ((vx_L * vx_L) + (vy_L * vy_L));
  double en_R =
      P_R / (gamma - 1) + 0.5 * rho_R * ((vx_R * vx_R) + (vy_R * vy_R));

  double rho_star = 0.5 * (rho_L + rho_R);

  double momx_star = 0.5 * (rho_L * vx_L + rho_R * vx_R);
  double momy_star = 0.5 * (rho_L * vy_L + rho_R * vy_R);
  double en_star = 0.5 * (en_L + en_R);

  double P_star =
      (gamma - 1) *
      (en_star -
       0.5 * ((momx_star * momx_star) + (momy_star * momy_star)) / rho_star);

  double flux_Mass = momx_star;
  double flux_Momx = (momx_star * momx_star / rho_star) + P_star;
  double flux_Momy = (momy_star * momx_star / rho_star);
  double flux_Energy = momx_star / rho_star * (en_star + P_star);

  double C_L = sqrt(gamma * P_L / rho_L) + fabs(vx_L);
  double C_R = sqrt(gamma * P_R / rho_R) + fabs(vx_R);
  double C = C_L > C_R ? C_L : C_R;

  flux_Mass -= (C * 0.5 * (rho_L - rho_R));
  flux_Momx -= (C * 0.5 * (rho_L * vx_L - rho_R * vx_R));
  flux_Momy -= (C * 0.5 * (rho_L * vy_L - rho_R * vy_R));
  flux_Energy -= (C * 0.5 * (en_L - en_R));

  *flux_Mass_o = flux_Mass;
  *flux_Momx_o = flux_Momx;
  *flux_Momy_o = flux_Momy;
  *flux_Energy_o = flux_Energy;
}

void linspace(double start, double stop, size_t num, double *output) {
  double delta = (stop - start) / (num - 1);
  size_t count = 0;
  double curr = start;

  while (count < num) {
    output[count] = curr;
    curr += delta;
    count++;
  }
}

void meshgrid(double *xv, double *yv, size_t N, double *x_res, double *y_res) {
  for (size_t y = 0; y < N; ++y) {
    for (size_t x = 0; x < N; ++x) {
      x_res[y * N + x] = xv[x];
    }
  }

  for (size_t x = 0; x < N; ++x) {
    for (size_t y = 0; y < N; ++y) {
      y_res[y * N + x] = yv[y];
    }
  }
}

void init_cells(double *X, double *Y, double gamma, double vol, double w0,
                double sigma, size_t N, double *Mass_o, double *Momx_o,
                double *Momy_o, double *Energy_o) {
  for (size_t i = 0; i < N * N; ++i) {
    double val = fabs(Y[i] - 0.5) < 0.25 ? 1.0 : 0.0;

    double rho = val + 1;
    double vx = val - 0.5;
    double vy = w0 * sin(4 * M_PI * X[i]) *
                (exp(-1 * (pow(Y[i] - 0.25, 2) / (2 * pow(sigma, 2)))) +
                 exp(-1 * (pow(Y[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
    double P = 2.5;

    double mass, momx, momy, energy;

    getConserved(rho, vx, vy, P, gamma, vol, &mass, &momx, &momy, &energy);
    Mass_o[i] = mass;
    Momx_o[i] = momx;
    Momy_o[i] = momy;
    Energy_o[i] = energy;
  }
}

void printMatrix(double *M, int N) {
  for (int y = 0; y < N; ++y) {
    for (int x = 0; x < N; ++x) {
      printf("%lf ", M[y * N + x]);
    }
    printf("\n");
  }
  printf("\n");
}

void printM(double *M, int width, int height) {
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      printf("%lf ", M[y * width + x]);
    }
    printf("\n");
  }
  printf("\n");
}

void simloop(int n) {
  // Simulation parameters
}

void getBlock(double *arr, int arr_w, int x, int y, int xdim, int ydim,
              double *out) {
  int i = y * arr_w + x;
  int j = 0;

  while (j < xdim * ydim) {
    for (int k = 0; k < xdim; ++k) {
      out[j++] = arr[i++];
    }
    i += arr_w - xdim;
  }
}

void sendRecieveGhostAllDir(double *sendbuf, int blocksize, int up, int down,
                            int left, int right, int blockWidth,
                            int blockHeight, double *up_o, double *down_o,
                            double *left_o, double *right_o) {
  int width = blockWidth + 2;
  int height = blockHeight + 2;

  double *temp_above = (double *)malloc(blockWidth * sizeof(double));
  double *temp_below = (double *)malloc(blockWidth * sizeof(double));
  double *temp_left = (double *)malloc(blockHeight * sizeof(double));
  double *temp_right = (double *)malloc(blockHeight * sizeof(double));

  double *left_temp = (double *)malloc(blockHeight * sizeof(double));
  double *right_temp = (double *)malloc(blockHeight * sizeof(double));
  getBlock(sendbuf, blockWidth + 2, 1, 1, 1, blockHeight, left_temp);
  getBlock(sendbuf, blockWidth + 2, blockWidth, 1, 1, blockHeight, right_temp);

  // sending ghost values
  MPI_Request reqs[8];
  MPI_Isend(&sendbuf[1 * (blockWidth + 2) + 1], blockWidth, MPI_DOUBLE, up, 1,
            MPI_COMM_WORLD, &reqs[0]);
  MPI_Isend(&sendbuf[blockHeight * (blockWidth + 2) + 1], blockWidth,
            MPI_DOUBLE, down, 2, MPI_COMM_WORLD, &reqs[1]);
  MPI_Isend(left_temp, blockHeight, MPI_DOUBLE, right, 1, MPI_COMM_WORLD,
            &reqs[2]);
  MPI_Isend(right_temp, blockHeight, MPI_DOUBLE, left, 2, MPI_COMM_WORLD,
            &reqs[3]);

  MPI_Request recvreqs[4];
  MPI_Irecv(temp_above, blockWidth, MPI_DOUBLE, down, 1, MPI_COMM_WORLD,
            &reqs[4]);
  MPI_Irecv(temp_below, blockWidth, MPI_DOUBLE, up, 2, MPI_COMM_WORLD,
            &reqs[5]);
  MPI_Irecv(temp_left, blockHeight, MPI_DOUBLE, right, 2, MPI_COMM_WORLD,
            &reqs[6]);
  MPI_Irecv(temp_right, blockHeight, MPI_DOUBLE, left, 1, MPI_COMM_WORLD,
            &reqs[7]);

  MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
}

void getGhostValsAllDir(double * f, int blockWidth, int blockHeight, int up, int down, int left, int right, double * above_o, double * below_o, double * left_o, double * right_o) {

    MPI_Request reqs[8];

    double * left_temp = (double *) malloc(sizeof(double) * blockHeight + 2);
    double * right_temp = (double *) malloc(sizeof(double) * blockHeight + 2);

    for(int i = 0; i < blockHeight; ++i) {
      left_temp[i + 1] = f[(1 + i) * (blockWidth + 2) + 1];
      right_temp[i + 1] = f[(1 + i) * (blockWidth + 2) + blockWidth];
    }

    left_temp[0] = left_temp[blockHeight + 1] = right_temp[0] = right_temp[blockHeight+ 1] = left_temp[1];

    MPI_Isend(&f[blockWidth + 2 + 1], blockWidth, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &reqs[0]);
    MPI_Isend(&f[(blockHeight) * (blockWidth + 2) + 1], blockWidth, MPI_DOUBLE, down, 2, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(left_temp, blockHeight + 2, MPI_DOUBLE, right, 3, MPI_COMM_WORLD, &reqs[2]);
    MPI_Isend(right_temp, blockHeight + 2, MPI_DOUBLE, left, 4, MPI_COMM_WORLD, &reqs[3]);

    MPI_Irecv(below_o, blockWidth, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, &reqs[5]);
    MPI_Irecv(above_o, blockWidth, MPI_DOUBLE, up, 2, MPI_COMM_WORLD, &reqs[4]);
    MPI_Irecv(right_o, blockHeight + 2, MPI_DOUBLE, left, 3, MPI_COMM_WORLD, &reqs[6]);
    MPI_Irecv(left_o, blockHeight + 2, MPI_DOUBLE, right, 4, MPI_COMM_WORLD, &reqs[7]);

    MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

    free(left_temp);
    free(right_temp);
}

void rollLeft(double * f, int blockWidth, int blockHeight, int blocksize, int rank_left, int rank_right) {
    double * col_temp = (double *) malloc(sizeof(double) * blockHeight);
    double * replace_col= (double *) malloc(sizeof(double) * blockHeight);
    for(int i = 0; i < blockHeight; ++i) { col_temp[i] = f[blockWidth * i]; }

    MPI_Request req;
    MPI_Isend(col_temp, blockHeight, MPI_DOUBLE, rank_left, 1, MPI_COMM_WORLD, &req);
    MPI_Recv(replace_col, blockWidth, MPI_DOUBLE, rank_right, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Wait(&req, MPI_STATUS_IGNORE);

    for(int y = 0; y < blockHeight; ++y) {
      for(int x = 1; x < blockWidth; ++x) {
        f[y * blockWidth + x - 1] = f[y * blockWidth + x];
      }
    }

    for(int i = 0; i < blockHeight; ++i) { f[i * blockWidth + (blockWidth - 1)] = replace_col[i]; }

    // roll matrix up
    free(col_temp);
    free(replace_col);
}

void rollRight(double * f, int blockWidth, int blockHeight, int blocksize, int rank_left, int rank_right) {
    
    double * col_temp = (double *) malloc(sizeof(double) * blockHeight);
    double * replace_col= (double *) malloc(sizeof(double) * blockHeight);

    for(int i = 0; i < blockHeight; ++i) { col_temp[i] = f[blockWidth * i + blockWidth - 1]; }
    

    MPI_Request req;
    MPI_Isend(col_temp, blockHeight, MPI_DOUBLE, rank_right, 1, MPI_COMM_WORLD, &req);
    MPI_Recv(replace_col, blockWidth, MPI_DOUBLE, rank_left, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Wait(&req, MPI_STATUS_IGNORE);

    // printf("after wait\n");

    for(int y = 0; y < blockHeight; ++y) {
      for(int x = blockWidth - 1; x > 0; --x) {
        f[x] = f[y * blockWidth + x - 1];
      }
    }

    for(int i = 0; i < blockHeight; ++i) { f[i * blockWidth] = replace_col[i]; }

    free(col_temp);
    free(replace_col);
}

void rollUp(double * f, int blockWidth, int blockHeight, int blocksize, int rank_up, int rank_down) {
  double * row_temp = (double *) malloc(sizeof(double) * blockWidth);

  MPI_Request req;

  MPI_Isend(f, blockWidth, MPI_DOUBLE, rank_up, 1, MPI_COMM_WORLD, &req);
  MPI_Recv(row_temp, blockWidth, MPI_DOUBLE, rank_down, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);

  // roll matrix up
  memcpy(f, &f[blockWidth], (blocksize - blockWidth) * sizeof(double));
  memcpy(&f[(blockHeight - 1) * blockWidth], row_temp, sizeof(double) * blockWidth);
  free(row_temp);
}

void rollDown(double * f, int blockWidth, int blockHeight, int blocksize, int rank_down, int rank_up) {
  double * row_temp = (double *) malloc(sizeof(double) * blockWidth);

  MPI_Request req;

  MPI_Isend(&f[(blockHeight - 1) * blockWidth], blockWidth, MPI_DOUBLE, rank_down, 1, MPI_COMM_WORLD, &req);
  MPI_Recv(row_temp, blockWidth, MPI_DOUBLE, rank_up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&req, MPI_STATUS_IGNORE);

  // roll matrix down
  memcpy(&f[blockWidth],f, (blocksize - blockWidth) * sizeof(double));
  memcpy(f, row_temp, sizeof(double) * blockWidth);
  free(row_temp);
}

int main(int argc, char *argv[]) {
  int rank, size, i, provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Initial conditions
  size_t N = 4;
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0;
  double courant_fac = 0.4;
  double t = 0.0;
  double tEnd = 2.0;
  double tOut = 0.5;
  double dx = boxsize / (double)N;
  double vol = dx * dx;
  double w0 = 0.1;
  double sigma = 0.05 / sqrt(2.0);

  // Only root needs a buffer to accumulate
  // values of rho, i.e result of the computation
  double *global_rho;
  if (rank == 0) {
    global_rho = (double *)malloc(N * N * sizeof(double));
  }

  // Check that work can be divided properly
  if ((N * N) % size != 0) {
    printf("N*N needs to be divisible by size\n");
    return 0;
  }

  // how many cells each worker is responsible for
  int blocksize = N * N / size;

  // height and width of each block in for all processes
  int blockWidth = blocksize > N ? N : blocksize;
  int blockHeight = blocksize < N ? 1 : blocksize / N;

  int nw = N / blockWidth;
  int nh = N / blockHeight;

  int l_y = rank / (N / blockWidth);
  int l_x = (rank) % (N / blockWidth);

  int rank_up = ((l_y + nh - 1) % nh) * (N / blockWidth) + l_x;
  int rank_down = ((l_y + 1) % nh) * (N / blockWidth) + l_x;

  int rank_left = l_y * (N / blockWidth) + ((l_x + 1) % nw);
  int rank_right = l_y * (N / blockWidth) + ((l_x + nw - 1) % nw);

  // Allocate memory for each process to store their unique values
  double *l_Mass = (double *)malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));
  double *l_Momx = (double *)malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));
  double *l_Momy = (double *)malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));
  double *l_Energy = (double *)malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));

  // // space for local values including ghost cells
  double *l_rho = (double *) malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));
  double *l_vx = (double *) malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));
  double *l_vy = (double *) malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));
  double *l_P = (double *) malloc(sizeof(double) * (blockWidth + 2) * (blockHeight + 2));

  double *l_rhodx = (double *) malloc(sizeof(double) * blocksize);
  double *l_rhody= (double *) malloc(sizeof(double) * blocksize);
  double *l_vxdx = (double *) malloc(sizeof(double) * blocksize);
  double *l_vxdy= (double *) malloc(sizeof(double) * blocksize);
  double *l_vydx = (double *) malloc(sizeof(double) * blocksize);
  double *l_vydy = (double *) malloc(sizeof(double) * blocksize);
  double *l_Pdx = (double *) malloc(sizeof(double) * blocksize);
  double *l_Pdy= (double *) malloc(sizeof(double) * blocksize);

  double *l_rhoprime= (double *) malloc(sizeof(double) * blocksize);
  double *l_vxprime = (double *) malloc(sizeof(double) * blocksize);
  double *l_vyprime = (double *) malloc(sizeof(double) * blocksize);
  double *l_Pprime = (double *) malloc(sizeof(double) * blocksize);

  double * rho_XL = (double *) malloc(sizeof(double) * blocksize);
  double * rho_XR = (double *) malloc(sizeof(double) * blocksize);
  double * rho_YL = (double *) malloc(sizeof(double) * blocksize);
  double * rho_YR = (double *) malloc(sizeof(double) * blocksize);

  double * vx_XL = (double *) malloc(sizeof(double) * blocksize);
  double * vx_XR = (double *) malloc(sizeof(double) * blocksize);
  double * vx_YL = (double *) malloc(sizeof(double) * blocksize);
  double * vx_YR = (double *) malloc(sizeof(double) * blocksize);

  double * vy_XL = (double *) malloc(sizeof(double) * blocksize);
  double * vy_XR = (double *) malloc(sizeof(double) * blocksize);
  double * vy_YL = (double *) malloc(sizeof(double) * blocksize);
  double * vy_YR = (double *) malloc(sizeof(double) * blocksize);

  double * P_XL = (double *) malloc(sizeof(double) * blocksize);
  double * P_XR = (double *) malloc(sizeof(double) * blocksize);
  double * P_YL = (double *) malloc(sizeof(double) * blocksize);
  double * P_YR = (double *) malloc(sizeof(double) * blocksize);

  double * flux_Mass_X = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Momx_X = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Momy_X = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Energy_X = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Mass_Y = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Momx_Y = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Momy_Y = (double *) malloc(sizeof(double) * blocksize);
  double * flux_Energy_Y = (double *) malloc(sizeof(double) * blocksize);


  // Generate and distribute inital vals for each process
  if (rank == 0) {
    // Generate initial values
    double *xlin = (double *)malloc(N * sizeof(double));
    double *X = (double *)malloc(N * N * sizeof(double));
    double *Y = (double *)malloc(N * N * sizeof(double));
    double *Mass = (double *)malloc(N * N * sizeof(double));
    double *Momx = (double *)malloc(N * N * sizeof(double));
    double *Momy = (double *)malloc(N * N * sizeof(double));
    double *Energy = (double *)malloc(N * N * sizeof(double));

    linspace(0.5 * dx, boxsize - 0.5 * dx, N, xlin);
    meshgrid(xlin, xlin, N, Y, X);

    init_cells(X, Y, gamma, vol, w0, sigma, N, Mass, Momx, Momy, Energy);

    MPI_Request req;
    // send initial cell values to all workers
    for (int i = 0; i < size; ++i) {
      // dont send cell values to ourselves
      if (i == 0) {
        memcpy(&l_Mass[blockWidth + 2 + 1], Mass, blocksize * sizeof(double));
        memcpy(&l_Momx[blockWidth + 2 + 1], Momx, blocksize * sizeof(double));
        memcpy(&l_Momy[blockWidth + 2 + 1], Momy, blocksize * sizeof(double));
        memcpy(&l_Energy[blockWidth + 2 + 1], Energy, blocksize * sizeof(double));
        continue;
      }

      MPI_Isend(Mass + (i * blocksize), blocksize, MPI_DOUBLE, i, 0,
                MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);

      MPI_Isend(Momx + (i * blocksize), blocksize, MPI_DOUBLE, i, 1,
                MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);

      MPI_Isend(Momy + (i * blocksize), blocksize, MPI_DOUBLE, i, 2,
                MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);

      MPI_Isend(Energy + (i * blocksize), blocksize, MPI_DOUBLE, i, 3,
                MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
    }

    free(xlin);
    free(X);
    free(Y);
    free(Mass);
    free(Momx);
    free(Momy);
    free(Energy);
  } else {
    MPI_Request requests[4];

    MPI_Irecv(&l_Mass[blockWidth + 2 + 1], blocksize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
              &requests[0]);
    MPI_Irecv(&l_Momx[blockWidth + 2 + 1], blocksize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,
              &requests[1]);
    MPI_Irecv(&l_Momy[blockWidth + 2 + 1], blocksize, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,
              &requests[2]);
    MPI_Irecv(&l_Energy[blockWidth + 2 + 1], blocksize, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD,
              &requests[3]);

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
  }

  size_t outputCount = 1;

  // Cellwise loop for multithreading, swap for and while for multi-threading
  while (t < tEnd) {
    printf("%lf / %lf\n", t, tEnd);
    double dt = 999999999999999999.0;

    // Send and recieve ghost vals for mass, momx, momy and energy
    double * left_temp = (double*)malloc(sizeof(double) * (blockHeight + 2)); 
    double * right_temp = (double*)malloc(sizeof(double) * (blockHeight + 2)); 

    getGhostValsAllDir(l_Mass, blockWidth, blockHeight, rank_up, rank_down, rank_left, rank_right, &l_Mass[1], &l_Mass[(blockHeight + 1) * (blockWidth + 2) + 1], left_temp, right_temp);
    for(int i = 0; i < blockHeight + 2; ++i) {
      l_Mass[(i) * (blockWidth + 2)] = left_temp[i];
      l_Mass[(i) * (blockWidth + 2) + blockWidth + 1] = right_temp[i];
    }

    getGhostValsAllDir(l_Momx, blockWidth, blockHeight, rank_up, rank_down, rank_left, rank_right, &l_Momx[1], &l_Momx[(blockHeight + 1) * (blockWidth + 2) + 1], left_temp, right_temp);
    for(int i = 0; i < blockHeight + 2; ++i) {
      l_Momx[(i) * (blockWidth + 2)] = left_temp[i];
      l_Momx[(i) * (blockWidth + 2) + blockWidth + 1] = right_temp[i];
    }

    getGhostValsAllDir(l_Momy, blockWidth, blockHeight, rank_up, rank_down, rank_left, rank_right, &l_Momy[1], &l_Momy[(blockHeight + 1) * (blockWidth + 2) + 1], left_temp, right_temp);
    for(int i = 0; i < blockHeight + 2; ++i) {
      l_Momy[(i) * (blockWidth + 2)] = left_temp[i];
      l_Momy[(i) * (blockWidth + 2) + blockWidth + 1] = right_temp[i];
    }

    getGhostValsAllDir(l_Energy, blockWidth, blockHeight, rank_up, rank_down, rank_left, rank_right, &l_Energy[1], &l_Energy[(blockHeight + 1) * (blockWidth + 2) + 1], left_temp, right_temp);
    for(int i = 0; i < blockHeight + 2; ++i) {
      l_Energy[(i) * (blockWidth + 2)] = left_temp[i];
      l_Energy[(i) * (blockWidth + 2) + blockWidth + 1] = right_temp[i];
    }

    free(left_temp);
    free(right_temp);

    for(int i = 0; i < (blockHeight + 2) * (blockWidth + 2); ++i) {
      getPrimitive(l_Mass[i], l_Momx[i], l_Momy[i], l_Energy[i], gamma, vol,
                       &l_rho[i], &l_vx[i], &l_vy[i], &l_P[i]);

        double dt_i =
            courant_fac * (dx / (sqrt(gamma * l_P[i] / l_rho[i]) +
                                 sqrt(pow(l_vx[i], 2) + pow(l_vy[i], 2))));  

        dt = dt < dt_i ? dt : dt_i;
    }
    // printf("rank %d calced primitives\n", rank);

    // sync dt value across workers
    if (rank == 0) {
      MPI_Request requests[size - 1];
      double dts[size - 1];

      // get local min dt values from other processes
      for (int r = 1; r < size; ++r) {
        MPI_Irecv(&dts[r - 1], 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD,
                  &requests[r - 1]);
      }

      MPI_Waitall(size - 1, requests, MPI_STATUSES_IGNORE);

      for (int i = 0; i < size - 1; ++i) {
        dt = dt < dts[i] ? dt : dts[i];
      }

    } else {
      MPI_Request req;
      MPI_Isend(&dt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &req);
      MPI_Request_free(&req);
    }

    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    bool plotThisTurn = false;
    if (t + dt > outputCount * tOut) {
      dt = outputCount * tOut - t;
      plotThisTurn = true;
    }
    t += dt;


    // If we should plot this turn all workers need to send their
    // local results to root so root can print a global resut
    if (1 /*plotThisTurn || t >= tEnd*/) {
      if (rank == 0) {
        memcpy(global_rho, &l_rho[blockWidth + 2 + 1], blocksize * sizeof(double));

        MPI_Request requests[size - 1];

        for (int i = 1; i < size; ++i) {
          MPI_Irecv(global_rho + (i * blocksize), blocksize, MPI_DOUBLE, i, 0,
                    MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(size - 1, requests, MPI_STATUSES_IGNORE);
        // printf("combined rho:\n");
        // printMatrix(global_rho, N);

      } else {
        MPI_Send(&l_rho[blockWidth + 2 + 1], blocksize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }
    }


    for(int y = 1; y < blockHeight + 1; ++y) {
      for(int x = 1; x < blockWidth + 1; ++x) {
        int i = y * (blockWidth + 2) + x - (blockWidth + 2 + 1);
        int j = y * (blockWidth + 2) + x;
        int up = (y - 1) * (blockWidth +2) + x;
        int down = (y + 1) * (blockWidth +2) + x;
        int left = y * (blockWidth +2) + x - 1;
        int right = y * (blockWidth +2) + x + 1;

        getGradient(l_rho[up], l_rho[down], l_rho[left], l_rho[right], dx, &l_rhodx[i], &l_rhody[i]);
        getGradient(l_vx[up], l_vx[down], l_vx[left], l_vx[right], dx, &l_vxdx[i], &l_vxdy[i]);
        getGradient(l_vy[up], l_vy[down], l_vy[left], l_vy[right], dx, &l_vydx[i], &l_vydy[i]);
        getGradient(l_P[up], l_P[down], l_P[left], l_P[right], dx, &l_Pdx[i], &l_Pdy[i]);
        
        if(rank == 0) 
        printf("%.10f \n", l_vydy[i]);

        l_rhoprime[i] = l_rho[j] - 0.5 * dt * ((l_vx[j] * l_rhodx[i]) + (l_rho[j] * l_vxdx[i]) + (l_vy[j] * l_rhody[i]) + (l_rho[j] * l_vydy[i]));
        l_vxprime[i] = l_vx[j] - (0.5 * dt) * (l_vx[j] * l_vxdx[i] + l_vy[j] * l_vxdy[i] + (1 / l_rho[j]) * l_Pdx[i]);
        l_vyprime[i] = l_vy[j] - (0.5 * dt) * (l_vx[j] * l_vydx[i] + l_vy[j] * l_vydy[i] + (1 / l_rho[j]) * l_Pdy[i]);
        l_Pprime[i] = l_P[j] - (0.5 * dt) * (gamma * l_P[j] * (l_vxdx[i] + l_vydy[i]) + l_vx[i] * l_Pdx[i] + l_vy[j] * l_Pdy[i]);
      }
    }
    sleep(rank);
    printf("rank %d has rhoprime: \n", rank);
    printM(l_rhoprime, blockWidth, blockHeight);
    sleep(10);


    for(int i = 0; i < blocksize; ++i) {
      extrapolateInSpaceToFace(l_rhoprime[i], l_rhodx[i], l_rhody[i], dx, &rho_XL[i], &rho_XR[i], &rho_YL[i], &rho_YR[i]);
      extrapolateInSpaceToFace(l_vxprime[i], l_vxdx[i], l_vxdy[i], dx, &vx_XL[i], &vx_XR[i], &vx_YL[i], &vx_YR[i]);
      extrapolateInSpaceToFace(l_vyprime[i], l_vydx[i], l_vydy[i], dx, &vy_XL[i], &vy_XR[i], &vy_YL[i], &vy_YR[i]);
      extrapolateInSpaceToFace(l_Pprime[i], l_Pdx[i], l_Pdy[i], dx, &P_XL[i], &P_XR[i], &P_YL[i], &P_YR[i]);
    }

    rollUp(rho_XL, blockWidth, blockHeight, blocksize, rank_up, rank_down);
    rollUp(vx_XL, blockWidth, blockHeight, blocksize, rank_up, rank_down);
    rollUp(vy_XL, blockWidth, blockHeight, blocksize, rank_up, rank_down);
    rollUp(P_XL, blockWidth, blockHeight, blocksize, rank_up, rank_down);
    

    printM(rho_XL, blockWidth, blockHeight);


    rollLeft(rho_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);
    rollLeft(vx_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);
    rollLeft(vy_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);
    rollLeft(P_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);

    for(int i = 0; i< blocksize; ++i) {
      getFlux(rho_XL[i], rho_XR[i], vx_XL[i], vx_XR[i], vy_XL[i], vy_XR[i], P_XL[i], P_XR[i], gamma, &flux_Mass_X[i], &flux_Momx_X[i], &flux_Momy_X[i], &flux_Energy_X[i]);
      getFlux(rho_YL[i], rho_YR[i], vy_YL[i], vy_YR[i], vx_YL[i], vx_YR[i], P_YL[i], P_YR[i], gamma, &flux_Mass_Y[i], &flux_Momy_Y[i], &flux_Momx_Y[i], &flux_Energy_Y[i]);
    }

    // get ghost vals for _Y from left
    // get ghost vals for _X from above

    double * flux_YL = (double *) malloc(sizeof(double) * blocksize);
    double * flux_XU = (double *) malloc(sizeof(double) * blocksize);

    memcpy(flux_YL, flux_Mass_Y, sizeof(double) * blocksize);
    memcpy(flux_XU, flux_Mass_X, sizeof(double) * blocksize);
    rollDown(flux_XU ,blockWidth, blockHeight, blocksize,rank_down, rank_up);
    rollRight(flux_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);

    for(int i = 0; i < blocksize; ++i) {
      l_Mass[i] = applyFluxes(l_Mass[i], flux_Mass_X[i], flux_XU[i], flux_YL[i], flux_Mass_Y[i],dx, dt);
    }

    memcpy(flux_YL, flux_Momx_Y, sizeof(double) * blocksize);
    memcpy(flux_XU, flux_Momx_X, sizeof(double) * blocksize);
    rollRight(flux_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);
    rollDown(flux_XU ,blockWidth, blockHeight, blocksize,rank_down, rank_up);

    for(int i = 0; i < blocksize; ++i) {
      l_Momx[i] = applyFluxes(l_Momx[i], flux_Momx_X[i], flux_XU[i], flux_YL[i], flux_Momx_Y[i],dx, dt);
    }

    memcpy(flux_YL, flux_Momy_Y, sizeof(double) * blocksize);
    memcpy(flux_XU, flux_Momy_X, sizeof(double) * blocksize);

    rollRight(flux_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);
    rollDown(flux_XU ,blockWidth, blockHeight, blocksize,rank_down, rank_up);

    for(int i = 0; i < blocksize; ++i) {
      l_Momy[i] = applyFluxes(l_Momy[i], flux_Momy_X[i], flux_XU[i], flux_YL[i], flux_Momy_Y[i],dx, dt);
    }

    memcpy(flux_YL, flux_Energy_Y, sizeof(double) * blocksize);
    memcpy(flux_XU, flux_Energy_X, sizeof(double) * blocksize);
    rollRight(flux_YL, blockWidth, blockHeight, blocksize, rank_left, rank_right);
    rollDown(flux_XU ,blockWidth, blockHeight, blocksize,rank_down, rank_up);

    for(int i = 0; i < blocksize; ++i) {
      l_Energy[i] = applyFluxes(l_Energy[i], flux_Energy_X[i], flux_XU[i], flux_YL[i], flux_Energy_Y[i],dx, dt);
    }

    free(flux_YL);
    free(flux_XU);
      
    printf("rank %d finished iteration\n", rank);
    sleep(5);
    }

  free(l_Mass); 
  free(l_Momx); 
  free(l_Momy); 
  free(l_Energy); 

  // // space for local values including ghost cells
  free(l_rho); 
  free(l_vx); 
  free(l_vy); 
  free(l_P); 

  free(l_rhodx); 
  free(l_rhody);
  free(l_vxdx); 
  free(l_vxdy);
  free(l_vydx); 
  free(l_vydy); 
  free(l_Pdx); 
  free(l_Pdy);

  free(l_rhoprime);
  free(l_vxprime); 
  free(l_vyprime); 
  free(l_Pprime); 

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
