// #include "finitevolume.h"
#include <omp.h>
#include <unistd.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <format>
#include <iostream>
#include <vector>

// #include "matrix.h"
// np.roll(f,-1,axis=0) = rollUp
// np.roll(f,1,axis=0) = rollDown
// np.roll(rho,-1,axis=1) = rollLeft
// np.roll(rho,1,axis=1) = rollRight
//

struct matrix {
  double *data;
  size_t rows;
  size_t cols;

  matrix(int rows, int cols) : rows(rows), cols(cols) {
    data = new double[rows * cols];
  }

  ~matrix() { delete[] data; }
};

struct cell {
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
};

std::array<double, 4> getConserved(double rho, double vx, double vy, double P,
                                   double gamma, double vol) {
  double Mass = rho * vol;
  double Momx = Mass * vx;
  double Momy = Mass * vy;
  double Energy = (P / (gamma - 1) + 0.5 * rho * ((vx * vx) + (vy * vy))) * vol;

  return {Mass, Momx, Momy, Energy};
}

std::array<double, 4> getPrimitive(double Mass, double Momx, double Momy,
                                   double Energy, double gamma, double vol) {
  double rho = Mass / vol;
  double vx = Momx / rho / vol;
  double vy = Momy / rho / vol;
  double P = (Energy / vol - 0.5 * rho * ((vx * vx) + (vy * vy))) * (gamma - 1);

  return {rho, vx, vy, P};
}

// auto [vy_dx, vy_dy] =
//     getGradient(cells[i].vy, cells[up].vy, cells[down].vy, cells[left].vy,
//                 cells[right].vy, dx);
std::array<double, 2> getGradient(double f, double f_down, double f_up,
                                  double f_right, double f_left, double dx) {
  double f_dx = (f_up - f_down) / (2 * dx);
  double f_dy = (f_left - f_right) / (2 * dx);

  return {f_dx, f_dy};
}

// std::array<mtx, 2> slopeLimit(const mtx &f, double dx, mtx f_dx, mtx f_dy) {
//   f_dx *= Matrix::clamp(((f - f.rollDown()) / dx) / (f_dx.zeroCheck(1.0e-8)),
//   0, 1); f_dx *= Matrix::clamp(((-1) * (f - f.rollUp()) / dx) /
//   (f_dx.zeroCheck(1.0e-8)), 0, 1);

//   f_dy *= Matrix::clamp(((f - f.rollRight()) / dx) /
//   (f_dy.zeroCheck(1.0e-8)), 0, 1); f_dy *= Matrix::clamp(((-1) * (f -
//   f.rollLeft()) / dx) / (f_dy.zeroCheck(1.0e-8)), 0, 1);

//   return {f_dx, f_dy};
// }

std::array<double, 4> extrapolateInSpaceToFace(double f, double f_dx,
                                               double f_dy, double dx) {
  double f_XL = f - f_dx * (dx / 2.0);
  // f_XL = f_XL.rollUp();
  double f_XR = f + f_dx * (dx / 2.0);

  double f_YL = f - f_dy * (dx / 2.0);
  // f_YL = f_YL.rollLeft();
  double f_YR = f + f_dy * (dx / 2.0);

  return {f_XL, f_XR, f_YL, f_YR};
}

double applyFluxes(double F, double flux_F_x, double flux_F_x_up,
                   double flux_F_y_left, double flux_F_y, double dx,
                   double dt) {
  return F + (dt * dx) * (-flux_F_x + flux_F_x_up - flux_F_y + flux_F_y_left);
}

std::array<double, 4> getFlux(double rho_L, double rho_R, double vx_L,
                              double vx_R, double vy_L, double vy_R, double P_L,
                              double P_R, double gamma) {
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
  double C = std::max(C_L, C_R);

  flux_Mass -= (C * 0.5 * (rho_L - rho_R));
  flux_Momx -= (C * 0.5 * (rho_L * vx_L - rho_R * vx_R));
  flux_Momy -= (C * 0.5 * (rho_L * vy_L - rho_R * vy_R));
  flux_Energy -= (C * 0.5 * (en_L - en_R));

  return {flux_Mass, flux_Momx, flux_Momy, flux_Energy};
}

std::vector<double> linspace(double start, double stop, size_t num) {
  std::vector<double> result;

  double delta = (stop - start) / (num - 1);
  double count = 0, curr = start;
  while (count < num) {
    result.push_back(curr);
    curr += delta;
    count++;
  }

  return result;
}

void meshgrid(const std::vector<double> &xv, const std::vector<double> &yv,
              matrix &x_res, matrix &y_res) {
  for (size_t y = 0; y < x_res.rows; ++y) {
    for (size_t x = 0; x < x_res.cols; ++x) {
      x_res.data[y * x_res.cols + x] = xv[x];
    }
  }

  for (size_t x = 0; x < y_res.cols; ++x) {
    for (size_t y = 0; y < y_res.rows; ++y) {
      y_res.data[y * y_res.cols + x] = yv[y];
    }
  }
}

void init_cells(std::vector<cell> &cells, const matrix &X, const matrix &Y,
                double gamma, double vol, double w0, double sigma) {
  for (size_t i = 0; i < cells.size(); ++i) {
    double val = fabs(Y.data[i] - 0.5) < 0.25 ? 1.0 : 0.0;

    double rho = val + 1;
    double vx = val - 0.5;
    double vy = w0 * sin(4 * M_PI * X.data[i]) *
                (exp(-1 * (pow(Y.data[i] - 0.25, 2) / (2 * pow(sigma, 2)))) +
                 exp(-1 * (pow(Y.data[i] - 0.75, 2) / (2 * pow(sigma, 2)))));
    double P = 2.5;

    auto [Mass, Momx, Momy, Energy] = getConserved(rho, vx, vy, P, gamma, vol);
    cells[i] = cell{Mass, Momx, Momy, Energy, gamma, vol};
  }
}

void printRho(std::vector<cell> &cells) {
  for (int i = 0; i < cells.size(); ++i) {
    printf("%lf ", cells[i].rho);
  }
  printf("\n");
}

void printMass(std::vector<cell> &cells) {
  for (int i = 0; i < cells.size(); ++i) {
    printf("%lf ", cells[i].Mass);
  }
  printf("\n");
}

void printRhodx(std::vector<cell> &cells, size_t N) {
  for (size_t y = 0; y < N; ++y) {
    for (size_t x = 0; x < N; ++x) {
      printf("%lf ", cells[y * N + x].rho_dx);
    }
    printf("\n");
  }
  printf("\n");
}

void printRhody(std::vector<cell> &cells, size_t N) {
  for (size_t y = 0; y < N; ++y) {
    for (size_t x = 0; x < N; ++x) {
      printf("%lf ", cells[y * N + x].rho_dy);
    }
    printf("\n");
  }
  printf("\n");
}

void printCell(cell &c) {
  printf("Mass %lf\n", c.Mass);
  printf("Momx %lf\n", c.Momx);
  printf("Momy %lf\n", c.Momy);
  printf("Energy %lf\n", c.Energy);
  printf("gamma %lf\n", c.gamma);
  printf("vol %lf\n", c.vol);

  printf("rho %lf\n", c.rho);
  printf("vx %lf\n", c.vx);
  // std::cout << std::format("{}", c.vy);
  printf("vy %lf\n", c.vy);
  printf("P); %lf\n", c.P);

  printf("rho_prime %lf\n", c.rho_prime);
  printf("vx_prime %lf\n", c.vx_prime);
  printf("vy_prime %lf\n", c.vy_prime);
  printf("P_prime %lf\n", c.P_prime);

  printf("rho_dx %lf\n", c.rho_dx);
  printf("rho_dy %lf\n", c.rho_dy);

  printf("vx_dx %lf\n", c.vx_dx);
  printf("vx_dy %lf\n", c.vx_dy);

  printf("vy_dx %lf\n", c.vy_dx);
  printf("vy_dy %lf\n", c.vy_dy);

  printf("P_dx %lf\n", c.P_dx);
  printf("P_dy %lf\n", c.P_dy);

  printf("rho_XL %lf\n", c.rho_XL);
  printf("rho_XR %lf\n", c.rho_XR);
  printf("rho_YL %lf\n", c.rho_YL);
  printf("rho_YR %lf\n", c.rho_YR);

  printf("vx_XL %lf\n", c.vx_XL);
  printf("vx_XR %lf\n", c.vx_XR);
  printf("vx_YL %lf\n", c.vx_YL);
  printf("vx_YR %lf\n", c.vx_YR);

  printf("vy_XL %lf\n", c.vy_XL);
  printf("vy_XR %lf\n", c.vy_XR);
  printf("vy_YL %lf\n", c.vy_YL);
  printf("vy_YR %lf\n", c.vy_YR);

  printf("P_XL %lf\n", c.P_XL);
  printf("P_XR %lf\n", c.P_XR);
  printf("P_YL %lf\n", c.P_YL);
  printf("P_YR %lf\n", c.P_YR);

  printf("flux_Mass_X %lf\n", c.flux_Mass_X);
  printf("flux_Mass_Y %lf\n", c.flux_Mass_Y);

  printf("flux_Momx_X %lf\n", c.flux_Momx_X);
  printf("flux_Momx_Y %lf\n", c.flux_Momx_Y);

  printf("flux_Momy_X %lf\n", c.flux_Momy_X);
  printf("flux_Momy_Y %lf\n", c.flux_Momy_Y);

  printf("flux_Energy_X %lf\n", c.flux_Energy_X);
  printf("flux_Energy_Y %lf\n", c.flux_Energy_Y);
}

void simloop() {
  // FILE *stream = fopen("output.bin", "wb");
  // TODO: Write N as first line of output.bin

  // Simulation parameters
  size_t N = 256;
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
  std::vector<double> xlin = linspace(0.5 * dx, boxsize - 0.5 * dx, N);

  // Initialize matrices
  matrix X(xlin.size(), xlin.size());
  matrix Y(xlin.size(), xlin.size());
  meshgrid(xlin, xlin, Y, X);

  // Initialize cells
  std::vector<cell> cells(N * N);
  cells.resize(N * N);
  init_cells(cells, X, Y, gamma, vol, w0, sigma);

  size_t outputCount = 1;

  // Cellwise loop for multithreading, swap for and while for multi-threading
  while (t < tEnd) {
    // Sync minimum dt from all cells and get nearby cells
    double dt = std::numeric_limits<double>::max();
    // printf("%lf / %lf\n", t, tEnd);

#pragma omp parallel for
    for (size_t i = 0; i < cells.size(); ++i) {
      auto [rho, vx, vy, P] =
          getPrimitive(cells[i].Mass, cells[i].Momx, cells[i].Momy,
                       cells[i].Energy, gamma, vol);
      cells[i].rho = rho;
      cells[i].vx = vx;
      cells[i].vy = vy;
      cells[i].P = P;

      double dt_i =
          courant_fac *
          (dx / (sqrt(gamma * P / rho) + sqrt(pow(vx, 2) + pow(vy, 2))));
#pragma omp critical
      { dt = std::min(dt, dt_i); }
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
      // rho.writeToFile(stream);
      // if (t >= 1.95) {
      // for (size_t i = 0; i < cells.size(); ++i) {
      //   printf("%lf ", cells[i].rho);
      // }
      // printf("\n");
      // }
    }

// Calculate next step for all cells
#pragma omp parallel for
    for (size_t i = 0; i < cells.size(); ++i) {
      // auto [rho, vx, vy, P] = getPrimitive(cells[i].Mass, cells[i].Momx,
      // cells[i].Momy, cells[i].Energy, gamma, vol);

      // Indexing
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t down = (y + 1) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;
      size_t right = y * N + (x + 1) % N;

      // Requires nearby cells
      auto [rho_dx, rho_dy] =
          getGradient(cells[i].rho, cells[up].rho, cells[down].rho,
                      cells[left].rho, cells[right].rho, dx);
      cells[i].rho_dx = rho_dx;
      cells[i].rho_dy = rho_dy;

      auto [vx_dx, vx_dy] =
          getGradient(cells[i].vx, cells[up].vx, cells[down].vx, cells[left].vx,
                      cells[right].vx, dx);
      cells[i].vx_dx = vx_dx;
      cells[i].vx_dy = vx_dy;

      auto [vy_dx, vy_dy] =
          getGradient(cells[i].vy, cells[up].vy, cells[down].vy, cells[left].vy,
                      cells[right].vy, dx);

      cells[i].vy_dx = vy_dx;
      cells[i].vy_dy = vy_dy;

      auto [P_dx, P_dy] = getGradient(cells[i].P, cells[up].P, cells[down].P,
                                      cells[left].P, cells[right].P, dx);
      cells[i].P_dx = P_dx;
      cells[i].P_dy = P_dy;

      // Extrapolate half step in time
      cells[i].rho_prime =
          cells[i].rho - 0.5 * dt *
                             ((cells[i].vx * rho_dx) + (cells[i].rho * vx_dx) +
                              (cells[i].vy * rho_dy) + (cells[i].rho * vy_dy));

      cells[i].vx_prime =
          cells[i].vx -
          (0.5 * dt) * (cells[i].vx * vx_dx + cells[i].vy * vx_dy +
                        (1 / cells[i].rho) * P_dx);
      cells[i].vy_prime =
          cells[i].vy -
          (0.5 * dt) * (cells[i].vx * vy_dx + cells[i].vy * vy_dy +
                        (1 / cells[i].rho) * P_dy);
      cells[i].P_prime =
          cells[i].P - (0.5 * dt) * (gamma * cells[i].P * (vx_dx + vy_dy) +
                                     cells[i].vx * P_dx + cells[i].vy * P_dy);
    }

#pragma omp parallel for
    for (size_t i = 0; i < cells.size(); ++i) {
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;

      auto [rho_XL, rho_XR, rho_YL, rho_YR] = extrapolateInSpaceToFace(
          cells[i].rho_prime, cells[i].rho_dx, cells[i].rho_dy, dx);
      cells[up].rho_XL = rho_XL;
      cells[left].rho_YL = rho_YL;
      cells[i].rho_XR = rho_XR;
      cells[i].rho_YR = rho_YR;

      auto [vx_XL, vx_XR, vx_YL, vx_YR] = extrapolateInSpaceToFace(
          cells[i].vx_prime, cells[i].vx_dx, cells[i].vx_dy, dx);
      cells[up].vx_XL = vx_XL;
      cells[left].vx_YL = vx_YL;
      cells[i].vx_XR = vx_XR;
      cells[i].vx_YR = vx_YR;

      auto [vy_XL, vy_XR, vy_YL, vy_YR] = extrapolateInSpaceToFace(
          cells[i].vy_prime, cells[i].vy_dx, cells[i].vy_dy, dx);
      cells[up].vy_XL = vy_XL;
      cells[left].vy_YL = vy_YL;
      cells[i].vy_XR = vy_XR;
      cells[i].vy_YR = vy_YR;

      auto [P_XL, P_XR, P_YL, P_YR] = extrapolateInSpaceToFace(
          cells[i].P_prime, cells[i].P_dx, cells[i].P_dy, dx);
      cells[up].P_XL = P_XL;
      cells[left].P_YL = P_YL;
      cells[i].P_XR = P_XR;
      cells[i].P_YR = P_YR;
    }

#pragma omp parallel for
    for (size_t i = 0; i < cells.size(); ++i) {
      // Get fluxes
      auto [flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X] = getFlux(
          cells[i].rho_XL, cells[i].rho_XR, cells[i].vx_XL, cells[i].vx_XR,
          cells[i].vy_XL, cells[i].vy_XR, cells[i].P_XL, cells[i].P_XR, gamma);

      cells[i].flux_Mass_X = flux_Mass_X;
      cells[i].flux_Momx_X = flux_Momx_X;
      cells[i].flux_Momy_X = flux_Momy_X;
      cells[i].flux_Energy_X = flux_Energy_X;

      auto [flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y] = getFlux(
          cells[i].rho_YL, cells[i].rho_YR, cells[i].vy_YL, cells[i].vy_YR,
          cells[i].vx_YL, cells[i].vx_YR, cells[i].P_YL, cells[i].P_YR, gamma);

      cells[i].flux_Mass_Y = flux_Mass_Y;
      cells[i].flux_Momx_Y = flux_Momx_Y;
      cells[i].flux_Momy_Y = flux_Momy_Y;
      cells[i].flux_Energy_Y = flux_Energy_Y;
      // Rolls down and right
    }

#pragma omp parallel for
    for (size_t i = 0; i < cells.size(); ++i) {
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
}

int main() {
  omp_set_num_threads(2);
  printf("threads: %d\n", omp_get_max_threads());

  simloop();
  return 0;
}
