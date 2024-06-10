#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdio.h>

// 12 flop
void getConserved(double rho, double vx, double vy, double P, double gamma, double vol, double *Mass_o, double *Momx_o, double *Momy_o,
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

// 14 flop
void getPrimitive(double Mass, double Momx, double Momy, double Energy, double gamma, double vol, double *rho_o, double *vx_o, double *vy_o,
                  double *P_o) {
  double rho = Mass / vol;
  double vx = Momx / rho / vol;
  double vy = Momy / rho / vol;
  double P = (Energy / vol - 0.5 * rho * ((vx * vx) + (vy * vy))) * (gamma - 1);

  *rho_o = rho;
  *vx_o = vx;
  *vy_o = vy;
  *P_o = P;
}

// 6 flop
void getGradient(double f, double f_down, double f_up, double f_right, double f_left, double dx, double *f_dx_o, double *f_dy_o) {
  double f_dx = (f_up - f_down) / (2 * dx);
  double f_dy = (f_left - f_right) / (2 * dx);

  *f_dx_o = f_dx;
  *f_dy_o = f_dy;
}

// 12 flop
void extrapolateInSpaceToFace(double f, double f_dx, double f_dy, double dx, double *f_XL_o, double *f_XR_o, double *f_YL_o, double *f_YR_o) {
  double f_XL = f - f_dx * (dx / 2.0);
  double f_XR = f + f_dx * (dx / 2.0);

  double f_YL = f - f_dy * (dx / 2.0);
  double f_YR = f + f_dy * (dx / 2.0);

  *f_XL_o = f_XL;
  *f_XR_o = f_XR;
  *f_YR_o = f_YR;
  *f_YL_o = f_YL;
}

// 7 flop
double applyFluxes(double F, double flux_F_x, double flux_F_x_up, double flux_F_y_left, double flux_F_y, double dx, double dt) {
  return F + (dt * dx) * (-flux_F_x + flux_F_x_up - flux_F_y + flux_F_y_left);
}

// 71 flop
void getFlux(double rho_L, double rho_R, double vx_L, double vx_R, double vy_L, double vy_R, double P_L, double P_R, double gamma,
             double *flux_Mass_o, double *flux_Momx_o, double *flux_Momy_o, double *flux_Energy_o) {
  double en_L = P_L / (gamma - 1) + 0.5 * rho_L * ((vx_L * vx_L) + (vy_L * vy_L));
  double en_R = P_R / (gamma - 1) + 0.5 * rho_R * ((vx_R * vx_R) + (vy_R * vy_R));

  double rho_star = 0.5 * (rho_L + rho_R);

  double momx_star = 0.5 * (rho_L * vx_L + rho_R * vx_R);
  double momy_star = 0.5 * (rho_L * vy_L + rho_R * vy_R);
  double en_star = 0.5 * (en_L + en_R);

  double P_star = (gamma - 1) * (en_star - 0.5 * ((momx_star * momx_star) + (momy_star * momy_star)) / rho_star);

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


void printM(double *a, int w, int h) {
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      printf("%lf ", a[y * w + x]);
    }
    printf("\n");
  }
  printf("\n");
}

void writeToFile(FILE* stream, double * m, int w, int h) {
  fwrite(&w, sizeof(int), 1, stream); 
  fwrite(&h, sizeof(int), 1, stream); 
  fwrite(m, sizeof(double), w * h, stream);
}
