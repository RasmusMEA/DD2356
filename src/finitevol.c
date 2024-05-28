#include <omp.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>
#include <stdio.h>

typedef struct cells {
  double* Mass;
  double* Momx;
  double* Momy;
  double* Energy;
  double* gamma;
  double* vol;

  double* rho;
  double* vx;
  double* vy;
  double* P;

  double* rho_prime;
  double* vx_prime;
  double* vy_prime;
  double* P_prime;

  double* rho_dx;
  double* rho_dy;

  double* vx_dx;
  double* vx_dy;

  double* vy_dx;
  double* vy_dy;

  double* P_dx;
  double* P_dy;

  double* rho_XL;
  double* rho_XR;
  double* rho_YL;
  double* rho_YR;

  double* vx_XL;
  double* vx_XR;
  double* vx_YL;
  double* vx_YR;

  double* vy_XL;
  double* vy_XR;
  double* vy_YL;
  double* vy_YR;

  double* P_XL;
  double* P_XR;
  double* P_YL;
  double* P_YR;

  double* flux_Mass_X;
  double* flux_Mass_Y;

  double* flux_Momx_X;
  double* flux_Momx_Y;

  double* flux_Momy_X;
  double* flux_Momy_Y;

  double*  flux_Energy_X;
  double* flux_Energy_Y;

} cells;


void freeCells(cells c) {
  free(c.Mass);
 free(c.Momx);
 free(c.Momy);
 free(c.Energy);
 free(c.gamma);
 free(c.vol);

 free(c.rho);
 free(c.vx);
 free(c.vy);
 free(c.P);

 free(c.rho_prime);
 free(c.vx_prime);
 free(c.vy_prime);
 free(c.P_prime);

 free(c.rho_dx);
 free(c.rho_dy);

 free(c.vx_dx);
 free(c.vx_dy);

 free(c.vy_dx);
 free(c.vy_dy);

 free(c.P_dx);
 free(c.P_dy);

 free(c.rho_XL);
 free(c.rho_XR);
 free(c.rho_YL);
 free(c.rho_YR);

 free(c.vx_XL);
 free(c.vx_XR);
 free(c.vx_YL);
 free(c.vx_YR);

 free(c.vy_XL);
 free(c.vy_XR);
 free(c.vy_YL);
 free(c.vy_YR);

 free(c.P_XL);
 free(c.P_XR);
 free(c.P_YL);
 free(c.P_YR);

 free(c.flux_Mass_X);
 free(c.flux_Mass_Y);

 free(c.flux_Momx_X);
 free(c.flux_Momx_Y);

 free(c.flux_Momy_X);
 free(c.flux_Momy_Y);

 free(c. flux_Energy_X);
 free(c.flux_Energy_Y);
}


void initCells(cells * c, int N) {
  c->Mass = (double *)malloc(sizeof(double) * N * N);
  c->Momx = (double *)malloc(sizeof(double) * N * N);
  c->Momy = (double *)malloc(sizeof(double) * N * N);
  c->Energy = (double *)malloc(sizeof(double) * N * N);
  c->gamma = (double *)malloc(sizeof(double) * N * N);
  c->vol = (double *)malloc(sizeof(double) * N * N);

  c->rho = (double *)malloc(sizeof(double) * N * N);
  c->vx = (double *)malloc(sizeof(double) * N * N);
  c->vy = (double *)malloc(sizeof(double) * N * N);
  c->P = (double *)malloc(sizeof(double) * N * N);

  c->rho_prime = (double *)malloc(sizeof(double) * N * N);
  c->vx_prime = (double *)malloc(sizeof(double) * N * N);
  c->vy_prime = (double *)malloc(sizeof(double) * N * N);
  c->P_prime = (double *)malloc(sizeof(double) * N * N);

  c->rho_dx = (double *)malloc(sizeof(double) * N * N);
  c->rho_dy = (double *)malloc(sizeof(double) * N * N);

  c->vx_dx = (double *)malloc(sizeof(double) * N * N);
  c->vx_dy = (double *)malloc(sizeof(double) * N * N);

  c->vy_dx = (double *)malloc(sizeof(double) * N * N);
  c->vy_dy = (double *)malloc(sizeof(double) * N * N);

  c->P_dx = (double *)malloc(sizeof(double) * N * N);
  c->P_dy = (double *)malloc(sizeof(double) * N * N);

  c->rho_XL = (double *)malloc(sizeof(double) * N * N);
  c->rho_XR = (double *)malloc(sizeof(double) * N * N);
  c->rho_YL = (double *)malloc(sizeof(double) * N * N);
  c->rho_YR = (double *)malloc(sizeof(double) * N * N);

  c->vx_XL = (double *)malloc(sizeof(double) * N * N);
  c->vx_XR = (double *)malloc(sizeof(double) * N * N);
  c->vx_YL = (double *)malloc(sizeof(double) * N * N);
  c->vx_YR = (double *)malloc(sizeof(double) * N * N);

  c->vy_XL = (double *)malloc(sizeof(double) * N * N);
  c->vy_XR = (double *)malloc(sizeof(double) * N * N);
  c->vy_YL = (double *)malloc(sizeof(double) * N * N);
  c->vy_YR = (double *)malloc(sizeof(double) * N * N);

  c->P_XL = (double *)malloc(sizeof(double) * N * N);
  c->P_XR = (double *)malloc(sizeof(double) * N * N);
  c->P_YL = (double *)malloc(sizeof(double) * N * N);
  c->P_YR = (double *)malloc(sizeof(double) * N * N);

  c->flux_Mass_X = (double *)malloc(sizeof(double) * N * N);
  c->flux_Mass_Y = (double *)malloc(sizeof(double) * N * N);

  c->flux_Momx_X = (double *)malloc(sizeof(double) * N * N);
  c->flux_Momx_Y = (double *)malloc(sizeof(double) * N * N);

  c->flux_Momy_X = (double *)malloc(sizeof(double) * N * N);
  c->flux_Momy_Y = (double *)malloc(sizeof(double) * N * N);

   c->flux_Energy_X = (double *)malloc(sizeof(double) * N * N);
  c->flux_Energy_Y = (double *)malloc(sizeof(double) * N * N);
}

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

void getGradient(double f, double f_down, double f_up, double f_right,
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

void init_cells(cells *c, double *X, double *Y, double gamma, double vol,
                double w0, double sigma, size_t N) {
  for (size_t i = 0; i < N * N; ++i) {
    double val = fabs(Y[i] - 0.5) < 0.25 ? 1.0 : 0.0;

    double rho = val + 1;
    double vx = val - 0.5;
    double vy = w0 * sin(4 * M_PI * X[i]) *
                (exp(-1 * (pow(Y[i] - 0.25, 2) / (2 * pow(sigma, 2)))) +
                 exp(-1 * (pow(Y[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
    double P = 2.5;

    getConserved(rho, vx, vy, P, gamma, vol, &c->Mass[i], &c->Momx[i], &c->Momy[i],
                 &c->Energy[i]);

  }
}

void printRho(cells *c, int N) {
  for (int i = 0; i < N * N; ++i) {
    printf("%lf ", c->rho[i]);
  }
  printf("\n");
}

void simloop(int n) {
  // Simulation parameters
  size_t N = n;
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0;
  double courant_fac = 0.4;
  double t = 0.0;
  double tEnd = 2.0;
  double tOut = 0.02;
  bool useSlopeLimiting = false;
  double plotRealTime = true;

  // Initial conditions
  double w0 = 0.1;
  double sigma = 0.05 / sqrt(2.0);

  // Grid setup
  double dx = boxsize / (double)N;
  double vol = dx * dx;

  double *xlin = (double *)malloc(N * sizeof(double));
  double *X = (double *)malloc(N * N * sizeof(double));
  double *Y = (double *)malloc(N * N * sizeof(double));

  linspace(0.5 * dx, boxsize - 0.5 * dx, N, xlin);

  // Initialize matrices
  meshgrid(xlin, xlin, N, Y, X);

  // Initialize cells
  cells c; 
  initCells(&c, N);
  init_cells(&c, X, Y, gamma, vol, w0, sigma, N);

  free(xlin);
  free(X);
  free(Y);

  size_t outputCount = 1;



  // Cellwise loop for multithreading, swap for and while for multi-threading
  while (t < tEnd) {
    // Sync minimum dt from all cells and get nearby cells
    double dt = 999999999999999999.0;
    // printf("%lf / %lf\n", t, tEnd);

#pragma omp parallel
    {
      int arr_sec = N * N / omp_get_max_threads();
      int arr_extra = (N * N) % omp_get_max_threads();

      int n = omp_get_thread_num();
      int start = n * arr_sec + (n < arr_extra ? n : arr_extra);
      int end = start + arr_sec + (n < arr_extra ? 1 : 0);

      for (size_t i = start; i < end; ++i) {

        getPrimitive(c.Mass[i], c.Momx[i], c.Momy[i],
                     c.Energy[i], gamma, vol, &c.rho[i], &c.vx[i],
                     &c.vy[i], &c.P[i]);

        double dt_i = courant_fac *
                      (dx / (sqrt(gamma * c.P[i] / c.rho[i]) +
                             sqrt(pow(c.vx[i], 2) + pow(c.vy[i], 2))));

#pragma omp critical
        { dt = dt < dt_i ? dt : dt_i; }
      }

#pragma omp barrier

#pragma omp single
      {
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
          // if (t >= 1.95) {
          // printRho(cells, N);
          // }
        }
      }

      // Calculate next step for all cells
      for (size_t i = start; i < end; ++i) {
        // Indexing
        size_t y = i / N;
        size_t x = i % N;
        size_t up = (y - 1 + N) % N * N + x;
        size_t down = (y + 1) % N * N + x;
        size_t left = y * N + (x - 1 + N) % N;
        size_t right = y * N + (x + 1) % N;

        // Requires nearby cells
        getGradient(c.rho[i], c.rho[up], c.rho[down],
                    c.rho[left], c.rho[right], dx, &c.rho_dx[i],
                    &c.rho_dy[i]);

        getGradient(c.vx[i], c.vx[up], c.vx[down], c.vx[left],
                    c.vx[right], dx, &c.vx_dx[i], &c.vx_dy[i]);

        getGradient(c.vy[i], c.vy[up], c.vy[down], c.vy[left],
                    c.vy[right], dx, &c.vy_dx[i], &c.vy_dy[i]);

        getGradient(c.P[i], c.P[up], c.P[down], c.P[left],
                    c.P[right], dx, &c.P_dx[i], &c.P_dy[i]);

        // Extrapolate half step in time
        c.rho_prime[i] =
            c.rho[i] - 0.5 * dt *
                               ((c.vx[i] * c.rho_dx[i]) +
                                (c.rho[i] * c.vx_dx[i]) +
                                (c.vy[i] * c.rho_dy[i]) +
                                (c.rho[i] * c.vy_dy[i]));

        c.vx_prime[i] =
            c.vx[i] - (0.5 * dt) * (c.vx[i] * c.vx_dx[i] +
                                        c.vy[i] * c.vx_dy[i] +
                                        (1 / c.rho[i]) * c.P_dx[i]);
        c.vy_prime[i] =
            c.vy[i] - (0.5 * dt) * (c.vx[i] * c.vy_dx[i] +
                                        c.vy[i] * c.vy_dy[i] +
                                        (1 / c.rho[i]) * c.P_dy[i]);
        c.P_prime[i] =
            c.P[i] -
            (0.5 * dt) *
                (gamma * c.P[i] * (c.vx_dx[i] + c.vy_dy[i]) +
                 c.vx[i] * c.P_dx[i] + c.vy[i] * c.P_dy[i]);
      }

#pragma omp barrier

      for (size_t i = start; i < end; ++i) {
        size_t y = i / N;
        size_t x = i % N;
        size_t up = (y - 1 + N) % N * N + x;
        size_t left = y * N + (x - 1 + N) % N;

        extrapolateInSpaceToFace(c.rho_prime[i], c.rho_dx[i],
                                 c.rho_dy[i], dx, &c.rho_XL[up],
                                 &c.rho_XR[i], &c.rho_YL[left],
                                 &c.rho_YR[i]);

        extrapolateInSpaceToFace(c.vx_prime[i], c.vx_dx[i],
                                 c.vx_dy[i], dx, &c.vx_XL[up],
                                 &c.vx_XR[i], &c.vx_YL[left],
                                 &c.vx_YR[i]);

        extrapolateInSpaceToFace(c.vy_prime[i], c.vy_dx[i],
                                 c.vy_dy[i], dx, &c.vy_XL[up],
                                 &c.vy_XR[i], &c.vy_YL[left],
                                 &c.vy_YR[i]);

        extrapolateInSpaceToFace(c.P_prime[i], c.P_dx[i], c.P_dy[i],
                                 dx, &c.P_XL[up], &c.P_XR[i],
                                 &c.P_YL[left], &c.P_YR[i]);
      }

#pragma omp barrier
      for (size_t i = start; i < end; ++i) {
        getFlux(c.rho_XL[i], c.rho_XR[i], c.vx_XL[i],
                c.vx_XR[i], c.vy_XL[i], c.vy_XR[i], c.P_XL[i],
                c.P_XR[i], gamma, &c.flux_Mass_X[i],
                &c.flux_Momx_X[i], &c.flux_Momy_X[i],
                &c.flux_Energy_X[i]);

        getFlux(c.rho_YL[i], c.rho_YR[i], c.vy_YL[i],
                c.vy_YR[i], c.vx_YL[i], c.vx_YR[i], c.P_YL[i],
                c.P_YR[i], gamma, &c.flux_Mass_Y[i],
                &c.flux_Momy_Y[i], &c.flux_Momx_Y[i],
                &c.flux_Energy_Y[i]);
      }

#pragma omp barrier
      for (size_t i = start; i < end; ++i) {
        size_t y = i / N;
        size_t x = i % N;
        size_t up = (y - 1 + N) % N * N + x;
        size_t left = y * N + (x - 1 + N) % N;

        c.Mass[i] = applyFluxes(
            c.Mass[i], c.flux_Mass_X[i], c.flux_Mass_X[up],
            c.flux_Mass_Y[left], c.flux_Mass_Y[i], dx, dt);

        c.Momx[i] = applyFluxes(
            c.Momx[i], c.flux_Momx_X[i], c.flux_Momx_X[up],
            c.flux_Momx_Y[left], c.flux_Momx_Y[i], dx, dt);

        c.Momy[i] = applyFluxes(
            c.Momy[i], c.flux_Momy_X[i], c.flux_Momy_X[up],
            c.flux_Momy_Y[left], c.flux_Momy_Y[i], dx, dt);

        c.Energy[i] = applyFluxes(
            c.Energy[i], c.flux_Energy_X[i], c.flux_Energy_X[up],
            c.flux_Energy_Y[left], c.flux_Energy_Y[i], dx, dt);
      }
    }
  }

  // fclose(stream);
  printf("done!, outputcount: %zu\n", outputCount);
  freeCells(c);
}

int main(int argc, char* argv[]) {
  
  int a = 1;
  int N = 10;
  
  putenv( "OMP_WAIT_POLICY=ACTIVE" );

  if(argc < 3) {
    printf("Pass in thread number and resulotion as args\n");
    return 0;
  }
  else{
    a = atoi(argv[1]);
    N = atoi(argv[2]);
  }

  omp_set_num_threads(a);
  printf("threads: %d\n", omp_get_max_threads());
  

  double t1 = omp_get_wtime();
  simloop(N);
  double t2 = omp_get_wtime();
  printf("Time elapsed: %lf seconds\n", t2 - t1);
  return 0;
}
