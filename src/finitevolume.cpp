// #include "finitevolume.h"
#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>


// #include "matrix.h"
// np.roll(f,-1,axis=0) = rollUp
// np.roll(f,1,axis=0) = rollDown
// np.roll(rho,-1,axis=1) = rollLeft
// np.roll(rho,1,axis=1) = rollRight
//

struct matrix {
  double * data;
  size_t rows;
  size_t cols;

  matrix(int rows, int cols): rows(rows), cols(cols) {
    data = new double[rows * cols];
  }

  ~matrix() {
    delete [] data;
  }
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
};

std::array<double, 4> getConserved(double rho, double vx, double vy,
                                double P, double gamma, double vol) {
  double Mass = rho * vol;
  double Momx = Mass * vx;
  double Momy = Mass * vy;
  double Energy = (P / (gamma - 1) + 0.5 * rho * ((vx * vx) + (vy * vy))) * vol;

  return {Mass, Momx, Momy, Energy};
}

std::array<double, 4> getPrimitive(double Mass, double Momx,
                                double Momy, double Energy,
                                double gamma, double vol) {
  double rho = Mass / vol;
  double vx = Momx / rho / vol;
  double vy = Momy / rho / vol;
  double P = (Energy / vol - 0.5 * rho * ((vx * vx) + (vy * vy))) * (gamma - 1);

  return {rho, vx, vy, P};
}

std::array<double , 2> getGradient(double f, double f_up, double f_down, double f_left, double f_right, double dx) {
  double f_dx = (f_up - f_down )/ (2 * dx);
  double f_dy = (f_left - f_right) / (2 * dx);

  return {f_dx, f_dy};
}

// std::array<mtx, 2> slopeLimit(const mtx &f, double dx, mtx f_dx, mtx f_dy) {
//   f_dx *= Matrix::clamp(((f - f.rollDown()) / dx) / (f_dx.zeroCheck(1.0e-8)), 0, 1);
//   f_dx *= Matrix::clamp(((-1) * (f - f.rollUp()) / dx) / (f_dx.zeroCheck(1.0e-8)), 0, 1);

//   f_dy *= Matrix::clamp(((f - f.rollRight()) / dx) / (f_dy.zeroCheck(1.0e-8)), 0, 1);
//   f_dy *= Matrix::clamp(((-1) * (f - f.rollLeft()) / dx) / (f_dy.zeroCheck(1.0e-8)), 0, 1);

//   return {f_dx, f_dy};
// }

std::array<double, 4> extrapolateInSpaceToFace(double f, double f_dx, double f_dy, double dx) {
  double f_XL = f - f_dx * (dx / 2.0);
  double f_YL = f - f_dy * (dx / 2.0);

  // f_XL = f_XL.rollUp();
  // f_YL = f_YL.rollLeft();

  double f_XR = f + f_dx * (dx / 2.0);
  double f_YR = f + f_dy * (dx / 2.0);

  return {f_XL, f_XR, f_YL, f_YR};
}

double applyFluxes(double F, double flux_F_x, double flux_F_x_down, double flux_F_y_right, double flux_F_y, double dx, double dt) {
  return F + (dt * dx) * (- flux_F_x + flux_F_x_down - flux_F_y + flux_F_y_right);
}

std::array<double , 4> getFlux(double rho_L, double rho_R, double vx_L,
                              double vx_R, double vy_L, double vy_R,
                              double P_L, double P_R, double gamma) {
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

  double C_L = sqrt(gamma * P_L / rho_L) + abs(vx_L);
  double C_R = sqrt(gamma * P_R / rho_R) + abs(vx_R);
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

void meshgrid(const std::vector<double> &xv, const std::vector<double> &yv, matrix & x_res, matrix & y_res) {
  for (size_t y = 0; y < x_res.rows; ++y) {
    for (size_t x = 0; x < x_res.cols; ++x) {
      x_res.data[y * x_res.cols + x] = xv[x];
    }
  }

  for (size_t x = 0; x < y_res.cols; ++x) {
    for (size_t y = 0; y < y_res.rows; ++y) {
      y_res.data[y *  y_res.cols + x] = yv[y];
    }
  }
}

void init_cells(std::vector<cell> &cells, const matrix &X, const matrix &Y, double gamma, double vol, double w0, double sigma) {
  for(size_t i = 0; i < cells.size(); ++i) {
    double val = fabs(Y.data[i] - 0.5) < 0.25 ? 0.25 : (fabs(Y.data[i] - 0.5));
    
    double rho = val + 1;
    double vx = val - 0.5;
    double vy = w0 * sin(4 * M_PI * X.data[i]) * (
      exp(-1 * (pow(Y.data[i] - 0.25, 2)/(2 * pow(sigma, 2)))) + 
      exp(-1 * (pow(Y.data[i] - 0.75, 2)/(2 * pow(sigma, 2))))
    );
    double P = 2.5;
    
    auto [Mass, Momx, Momy, Energy] = getConserved(rho, vx, vy, P, gamma, vol);
    cells[i] = cell{Mass, Momx, Momy, Energy, gamma, vol};
  }
}

void simloop() {
  // FILE *stream = fopen("output.bin", "wb");
  // TODO: Write N as first line of output.bin

  // Simulation parameters
  size_t N = 8;
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
  meshgrid(xlin, xlin);
  
  // Initialize cells
  std::vector<cell> cells(N * N);
  cells.resize(N * N);
  init_cells(cells, X, Y, gamma, vol, w0, sigma);

  size_t outputCount = 1;

  // Cellwise loop for multithreading, swap for and while for multi-threading
  while (t < tEnd) {

    // Sync minimum dt from all cells and get nearby cells
    double dt = std::numeric_limits<double>::max();
    for (size_t i = 0; i < cells.size(); ++i) {
      printf("%lf / %lf\n", t, tEnd );
      auto [rho, vx, vy, P] = getPrimitive(cells[i].Mass, cells[i].Momx, cells[i].Momy, cells[i].Energy, gamma, vol);
      cells[i].rho = rho;
      cells[i].vx = vx;
      cells[i].vy = vy;
      cells[i].P = P;

      double dt_i = courant_fac * (dx / (sqrt(gamma * P / rho) + sqrt(pow(vx, 2) + pow(vy, 2))));
      dt = std::min(dt, dt_i);
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
      if(t >= 1.95) {
        for (size_t i = 0; i < cells.size(); ++i) {
          printf("%lf ", cells[i].rho);
        }
        printf("\n");
      }
    }

    // Calculate next step for all cells
    for (size_t i = 0; i < cells.size(); ++i) {
      // auto [rho, vx, vy, P] = getPrimitive(cells[i].Mass, cells[i].Momx, cells[i].Momy, cells[i].Energy, gamma, vol);

      // Indexing
      size_t y = i / N;
      size_t x = i % N;
      size_t up = (y - 1 + N) % N * N + x;
      size_t down = (y + 1) % N * N + x;
      size_t left = y * N + (x - 1 + N) % N;
      size_t right = y * N + (x + 1) % N;

      // Requires nearby cells
      auto [rho_dx, rho_dy] = getGradient(cells[i].rho, cells[up].rho, cells[down].rho, cells[left].rho, cells[right].rho, dx);
      auto [vx_dx, vx_dy] =  getGradient(cells[i].vx, cells[up].vx, cells[down].vx, cells[left].vx, cells[right].vx, dx);
      auto [vy_dx, vy_dy] = getGradient(cells[i].vy, cells[up].vy, cells[down].vy, cells[left].vy, cells[right].vy, dx);
      auto [P_dx, P_dy] = getGradient(cells[i].P, cells[up].P, cells[down].P, cells[left].P, cells[right].P, dx);

      // Extrapolate half step in time
      double rho_prime = cells[i].rho - (0.5 * dt) * (cells[i].vx * rho_dx + cells[i].rho * vx_dx + cells[i].vy * rho_dy + cells[i].rho * vy_dy);
      double vx_prime = cells[i].vx - (0.5 * dt) * (cells[i].vx * vx_dx + cells[i].vy * vx_dy + (1 / cells[i].rho) * P_dx);
      double vy_prime = cells[i].vy - (0.5 * dt) * (cells[i].vx * vy_dx + cells[i].vy * vy_dy + (1 / cells[i].rho) * P_dy);
      double P_prime = cells[i].P - (0.5 * dt) * (gamma * cells[i].P * (vx_dx + vy_dy) + cells[i].vx * P_dx + cells[i].vy * P_dy);
    }

    for (size_t i = 0; i < cells.size(); ++i) {
      auto [rho_XL, rho_XR, rho_YL, rho_YR] = extrapolateInSpaceToFace(rho_prime, rho_dx, rho_dy, dx);
      auto [vx_XL, vx_XR, vx_YL, vx_YR] = extrapolateInSpaceToFace(vx_prime, vx_dx, vx_dy, dx);
      auto [vy_XL, vy_XR, vy_YL, vy_YR] = extrapolateInSpaceToFace(vy_prime, vy_dx, vy_dy, dx);
      auto [P_XL, P_XR, P_YL, P_YR] = extrapolateInSpaceToFace(P_prime, P_dx, P_dy, dx);

      // Roll left and up here

      // Get fluxes
      auto [flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X] = getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, gamma);
      auto [flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y] = getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, gamma);

      // Rolls down and right

      // Apply fluxes
      cells[i].Mass = applyFluxes(cells[i].Mass, flux_Mass_X, flux_Mass_Y, dx, dt);
      cells[i].Momx = applyFluxes(cells[i].Momx, flux_Momx_X, flux_Momx_Y, dx, dt);
      cells[i].Momy = applyFluxes(cells[i].Momy, flux_Momy_X, flux_Momy_Y, dx, dt);
      cells[i].Energy = applyFluxes(cells[i].Energy, flux_Energy_X, flux_Energy_Y, dx, dt);
    }
  }

  // while (t < tEnd) {
  //   printf("%lf / %lf\n", t, tEnd );
  //   auto [rho, vx, vy, P] = getPrimitive(Mass, Momx, Momy, Energy, gamma, vol);

  //   // Requires syncronization
  //   double dt = courant_fac * Matrix::min(dx / (Matrix::sqrt(gamma * P / rho) + Matrix::sqrt((vx * vx) + (vy * vy))));

  //   bool plotThisTurn = false;

  //   if (t + dt > outputCount * tOut) {
  //     dt = outputCount * tOut - t;
  //     plotThisTurn = true;
  //   }

  //   // Requires nearby cells
  //   auto rho_d = getGradient(rho, dx);
  //   auto vx_d = getGradient(vx, dx);
  //   auto vy_d = getGradient(vy, dx);
  //   auto P_d = getGradient(P, dx);

  //   // if (useSlopeLimiting) {
  //   //   rho_d = slopeLimit(rho, dx, rho_d[0], rho_d[1]);
  //   //   vx_d = slopeLimit(vx, dx, vx_d[0], vx_d[1]);
  //   //   vy_d = slopeLimit(vy, dx, vy_d[0], vy_d[1]);
  //   //   P_d = slopeLimit(P, dx, P_d[0], P_d[1]);
  //   // }

  //   double rho_dx = std::move(rho_d[0]);
  //   double rho_dy = std::move(rho_d[1]);

  //   double vx_dx = std::move(vx_d[0]);
  //   double vx_dy = std::move(vx_d[1]);

  //   double vy_dx = std::move(vy_d[0]);
  //   double vy_dy = std::move(vy_d[1]);

  //   double P_dx = std::move(P_d[0]);
  //   double P_dy = std::move(P_d[1]);

  //   // Requires nearby cells
  //   // auto [rho_dx, rho_dy] = getGradient(rho, dx);
  //   // auto [vx_dx, vx_dy] = getGradient(vx, dx);
  //   // auto [vy_dx, vy_dy] = getGradient(vy, dx);
  //   // auto [P_dx, P_dy] = getGradient(P, dx);

  //   double rho_prime = rho - (0.5 * dt) * (vx * rho_dx + rho * vx_dx + vy * rho_dy + rho * vy_dy);
  //   double vx_prime = vx - (0.5 * dt) * (vx * vx_dx + vy * vx_dy + (1 / rho) * P_dx);
  //   double vy_prime = vy - (0.5 * dt) * (vx * vy_dx + vy * vy_dy + (1 / rho) * P_dy);
  //   double P_prime = P - (0.5 * dt) * (gamma * P * (vx_dx + vy_dy) + vx * P_dx + vy * P_dy);

  //   // Rolls up and left
  //   auto [rho_XL, rho_XR, rho_YL, rho_YR] = extrapolateInSpaceToFace(rho_prime, rho_dx, rho_dy, dx);
  //   auto [vx_XL, vx_XR, vx_YL, vx_YR] = extrapolateInSpaceToFace(vx_prime, vx_dx, vx_dy, dx);
  //   auto [vy_XL, vy_XR, vy_YL, vy_YR] = extrapolateInSpaceToFace(vy_prime, vy_dx, vy_dy, dx);
  //   auto [P_XL, P_XR, P_YL, P_YR] = extrapolateInSpaceToFace(P_prime, P_dx, P_dy, dx);

  //   auto [flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X] = getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, gamma);
  //   auto [flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y] = getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, gamma);

  //   Mass = applyFluxes(Mass, flux_Mass_X, flux_Mass_Y, dx, dt);
  //   Momx = applyFluxes(Momx, flux_Momx_X, flux_Momx_Y, dx, dt);
  //   Momy = applyFluxes(Momy, flux_Momy_X, flux_Momy_Y, dx, dt);
  //   Energy = applyFluxes(Energy, flux_Energy_X, flux_Energy_Y, dx, dt);

  //   t += dt;

  //   if (plotThisTurn || t >= tEnd) {
  //     outputCount += 1;
  //     // rho.writeToFile(stream);
  //     if(t >= 1.95)
  //     rho.print();
  //   }
  // }

  // fclose(stream);
  printf("done!, outputcount: %zu\n", outputCount);
}
