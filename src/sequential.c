#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct cell {
  double Mass;
  double Momx;
  double Momy;
  double Energy;
  double gamma;
  double vol;

  double rho;
  double vx;
  double vy;
  double P;

  double rho_prime;
  double vx_prime;
  double vy_prime;
  double P_prime;

  double rho_dx;
  double rho_dy;

  double vx_dx;
  double vx_dy;

  double vy_dx;
  double vy_dy;

  double P_dx;
  double P_dy;

  double rho_XL;
  double rho_XR;
  double rho_YL;
  double rho_YR;

  double vx_XL;
  double vx_XR;
  double vx_YL;
  double vx_YR;

  double vy_XL;
  double vy_XR;
  double vy_YL;
  double vy_YR;

  double P_XL;
  double P_XR;
  double P_YL;
  double P_YR;

  double flux_Mass_X;
  double flux_Mass_Y;

  double flux_Momx_X;
  double flux_Momx_Y;

  double flux_Momy_X;
  double flux_Momy_Y;

  double flux_Energy_X;
  double flux_Energy_Y;

} cell;

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

void init_cells(cell *cells, double *X, double *Y, double gamma, double vol,
                double w0, double sigma, size_t N) {
  for (size_t i = 0; i < N * N; ++i) {
    double val = fabs(Y[i] - 0.5) < 0.25 ? 1.0 : 0.0;

    double rho = val + 1;
    double vx = val - 0.5;
    double vy = w0 * sin(4 * M_PI * X[i]) *
                (exp(-1 * (pow(Y[i] - 0.25, 2) / (2 * pow(sigma, 2)))) +
                 exp(-1 * (pow(Y[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
    double P = 2.5;

    cell c;
    getConserved(rho, vx, vy, P, gamma, vol, &c.Mass, &c.Momx, &c.Momy,
                 &c.Energy);
    cells[i] = c;
  }
}

void printRho(cell *cells, int N) {
  for (int i = 0; i < N * N; ++i) {
    printf("%lf ", cells[i].rho);
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
  cell *cells = malloc(N * N * sizeof(cell));
  init_cells(cells, X, Y, gamma, vol, w0, sigma, N);

  free(xlin);
  free(X);
  free(Y);

  size_t outputCount = 1;

  // Cellwise loop for multithreading, swap for and while for multi-threading
  while (t < tEnd) {
    double dt = 999999999999999999.0;

    for (size_t i = 0; i < N * N; ++i) {
      getPrimitive(cells[i].Mass, cells[i].Momx, cells[i].Momy, cells[i].Energy,
                   gamma, vol, &cells[i].rho, &cells[i].vx, &cells[i].vy,
                   &cells[i].P);

      double dt_i = courant_fac *
                    (dx / (sqrt(gamma * cells[i].P / cells[i].rho) +
                           sqrt(pow(cells[i].vx, 2) + pow(cells[i].vy, 2))));

      dt = dt < dt_i ? dt : dt_i;
    }

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

    // Calculate next step for all cells
    for (size_t i = 0; i < N * N; ++i) {
      // Indexing
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t down = (y + 1) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;
      size_t right = y * N + (x + 1) % N;

      // Requires nearby cells
      getGradient(cells[i].rho, cells[up].rho, cells[down].rho, cells[left].rho,
                  cells[right].rho, dx, &cells[i].rho_dx, &cells[i].rho_dy);

      getGradient(cells[i].vx, cells[up].vx, cells[down].vx, cells[left].vx,
                  cells[right].vx, dx, &cells[i].vx_dx, &cells[i].vx_dy);

      getGradient(cells[i].vy, cells[up].vy, cells[down].vy, cells[left].vy,
                  cells[right].vy, dx, &cells[i].vy_dx, &cells[i].vy_dy);

      getGradient(cells[i].P, cells[up].P, cells[down].P, cells[left].P,
                  cells[right].P, dx, &cells[i].P_dx, &cells[i].P_dy);

      // Extrapolate half step in time
      cells[i].rho_prime = cells[i].rho - 0.5 * dt *
                                              ((cells[i].vx * cells[i].rho_dx) +
                                               (cells[i].rho * cells[i].vx_dx) +
                                               (cells[i].vy * cells[i].rho_dy) +
                                               (cells[i].rho * cells[i].vy_dy));

      cells[i].vx_prime =
          cells[i].vx - (0.5 * dt) * (cells[i].vx * cells[i].vx_dx +
                                      cells[i].vy * cells[i].vx_dy +
                                      (1 / cells[i].rho) * cells[i].P_dx);
      cells[i].vy_prime =
          cells[i].vy - (0.5 * dt) * (cells[i].vx * cells[i].vy_dx +
                                      cells[i].vy * cells[i].vy_dy +
                                      (1 / cells[i].rho) * cells[i].P_dy);
      cells[i].P_prime =
          cells[i].P -
          (0.5 * dt) *
              (gamma * cells[i].P * (cells[i].vx_dx + cells[i].vy_dy) +
               cells[i].vx * cells[i].P_dx + cells[i].vy * cells[i].P_dy);
    }

    for (size_t i = 0; i < N * N; ++i) {
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;

      extrapolateInSpaceToFace(cells[i].rho_prime, cells[i].rho_dx,
                               cells[i].rho_dy, dx, &cells[up].rho_XL,
                               &cells[i].rho_XR, &cells[left].rho_YL,
                               &cells[i].rho_YR);

      extrapolateInSpaceToFace(cells[i].vx_prime, cells[i].vx_dx,
                               cells[i].vx_dy, dx, &cells[up].vx_XL,
                               &cells[i].vx_XR, &cells[left].vx_YL,
                               &cells[i].vx_YR);

      extrapolateInSpaceToFace(cells[i].vy_prime, cells[i].vy_dx,
                               cells[i].vy_dy, dx, &cells[up].vy_XL,
                               &cells[i].vy_XR, &cells[left].vy_YL,
                               &cells[i].vy_YR);

      extrapolateInSpaceToFace(cells[i].P_prime, cells[i].P_dx, cells[i].P_dy,
                               dx, &cells[up].P_XL, &cells[i].P_XR,
                               &cells[left].P_YL, &cells[i].P_YR);
    }

    for (size_t i = 0; i < N * N; ++i) {
      getFlux(cells[i].rho_XL, cells[i].rho_XR, cells[i].vx_XL, cells[i].vx_XR,
              cells[i].vy_XL, cells[i].vy_XR, cells[i].P_XL, cells[i].P_XR,
              gamma, &cells[i].flux_Mass_X, &cells[i].flux_Momx_X,
              &cells[i].flux_Momy_X, &cells[i].flux_Energy_X);

      getFlux(cells[i].rho_YL, cells[i].rho_YR, cells[i].vy_YL, cells[i].vy_YR,
              cells[i].vx_YL, cells[i].vx_YR, cells[i].P_YL, cells[i].P_YR,
              gamma, &cells[i].flux_Mass_Y, &cells[i].flux_Momy_Y,
              &cells[i].flux_Momx_Y, &cells[i].flux_Energy_Y);
    }

    for (size_t i = 0; i < N * N; ++i) {
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;

      cells[i].Mass = applyFluxes(
          cells[i].Mass, cells[i].flux_Mass_X, cells[up].flux_Mass_X,
          cells[left].flux_Mass_Y, cells[i].flux_Mass_Y, dx, dt);

      cells[i].Momx = applyFluxes(
          cells[i].Momx, cells[i].flux_Momx_X, cells[up].flux_Momx_X,
          cells[left].flux_Momx_Y, cells[i].flux_Momx_Y, dx, dt);

      cells[i].Momy = applyFluxes(
          cells[i].Momy, cells[i].flux_Momy_X, cells[up].flux_Momy_X,
          cells[left].flux_Momy_Y, cells[i].flux_Momy_Y, dx, dt);

      cells[i].Energy = applyFluxes(
          cells[i].Energy, cells[i].flux_Energy_X, cells[up].flux_Energy_X,
          cells[left].flux_Energy_Y, cells[i].flux_Energy_Y, dx, dt);
    }
  }

  // fclose(stream);
  printf("done!, outputcount: %zu\n", outputCount);
  free(cells);
}

int main(int argc, char *argv[]) {
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
