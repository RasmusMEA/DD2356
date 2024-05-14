#include "finitevolume.h"

#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>

#include "matrix.h"

// np.roll(f,-1,axis=0) = rollUp
// np.roll(f,1,axis=0) = rollDown
// np.roll(rho,-1,axis=1) = rollLeft
// np.roll(rho,1,axis=1) = rollRight

std::array<mtx, 4> getConserved(const mtx &rho, const mtx &vx, const mtx &vy,
                                const mtx &P, double gamma, double vol) {
  mtx Mass = rho * vol;
  mtx Momx = rho * vx * vol;
  mtx Momy = rho * vy * vol;
  mtx Energy = (P / (gamma - 1) + 0.5 * rho * ((vx * vx) + (vy * vy))) * vol;

  return {Mass, Momx, Momy, Energy};
}

std::array<mtx, 4> getPrimitive(const mtx &Mass, const mtx &Momx,
                                const mtx &Momy, const mtx &Energy,
                                double gamma, double vol) {
  mtx rho = Mass / vol;
  mtx vx = Momx / rho / vol;
  mtx vy = Momy / rho / vol;
  mtx P = (Energy / vol - 0.5 * rho * ((vx * vx) + (vy * vy))) * (gamma - 1);

  return {rho, vx, vy, P};
}

std::array<mtx, 2> getGradient(const mtx &f, double dx) {
  mtx f_dx = (f.rollUp() - f.rollDown()) / (2 * dx);
  mtx f_dy = (f.rollLeft() - f.rollRight()) / (2 * dx);

  return {f_dx, f_dy};
}

std::array<mtx, 2> slopeLimit(const mtx &f, double dx, mtx f_dx, mtx f_dy) {
  f_dx = Matrix::maximum(
             Matrix::minimum(
                 ((f - f.rollDown()) / dx) / (f_dx.zeroCheck(1.0e-8)), 1),
             0) *
         f_dx;

  f_dx = Matrix::maximum(
             Matrix::minimum(
                 ((-1) * (f - f.rollUp()) / dx) / (f_dx.zeroCheck(1.0e-8)), 1),
             0) *
         f_dx;

  f_dy = Matrix::maximum(
             Matrix::minimum(
                 ((f - f.rollRight()) / dx) / (f_dy.zeroCheck(1.0e-8)), 1),
             0) *
         f_dy;

  f_dy = Matrix::maximum(Matrix::minimum(((-1) * (f - f.rollLeft()) / dx) /
                                             (f_dy.zeroCheck(1.0e-8)),
                                         1),
                         0) *
         f_dy;

  return {f_dx, f_dy};
}

std::array<mtx, 4> extrapolateInSpaceToFace(const mtx &f, const mtx &f_dx,
                                            const mtx &f_dy, double dx) {
  mtx f_XL = f - f_dx * dx / 2.0;
  f_XL = f_XL.rollUp();
  mtx f_XR = f + f_dx * dx / 2.0;

  mtx f_YL = f - f_dy * dx / 2.0;
  f_YL = f_YL.rollLeft();
  mtx f_YR = f + f_dy * dx / 2.0;

  return {f_XL, f_XR, f_YL, f_YR};
}

mtx applyFluxes(mtx F, const mtx &flux_F_x, const mtx &flux_F_y, double dx,
                double dt) {
  F = F - dt * dx * flux_F_x;
  F = F + dt * dx * flux_F_x.rollDown();
  F = F - dt * dx * flux_F_y;
  F = F + dt * dx * flux_F_y.rollRight();

  return F;
}

std::array<mtx, 4> getFlux(const mtx &rho_L, const mtx &rho_R, const mtx &vx_L,
                           const mtx &vx_R, const mtx &vy_L, const mtx &vy_R,
                           const mtx &P_L, const mtx &P_R, double gamma) {
  mtx en_L = P_L / (gamma - 1) + 0.5 * rho_L * ((vx_L * vx_L) + (vy_L * vy_L));
  mtx en_R = P_R / (gamma - 1) + 0.5 * rho_R * ((vx_R * vx_R) + (vy_R * vy_R));

  mtx rho_star = 0.5 * (rho_L + rho_R);
  mtx momx_star = 0.5 * (rho_L * vx_L + rho_R * vx_R);
  mtx momy_star = 0.5 * (rho_L * vy_L + rho_R * vy_R);
  mtx en_star = 0.5 * (en_L + en_R);

  mtx P_star =
      (gamma - 1) *
      (en_star -
       0.5 * ((momx_star * momx_star) + (momy_star * momy_star)) / rho_star);

  mtx flux_Mass = momx_star;
  mtx flux_Momx = ((momx_star * momx_star) / rho_star) + P_star;
  mtx flux_Momy = (momx_star * momy_star) / rho_star;
  mtx flux_Energy = (en_star + P_star) * momx_star / rho_star;

  mtx C_L = Matrix::sqrt(gamma * P_L / rho_L) + Matrix::abs(vx_L);
  mtx C_R = Matrix::sqrt(gamma * P_R / rho_R) + Matrix::abs(vx_R);
  mtx C = Matrix::maximum(C_L, C_R);

  flux_Mass = flux_Mass - (C * 0.5 * (rho_L - rho_R));
  flux_Momx = flux_Momx - (C * 0.5 * (rho_L * vx_L - rho_R * vx_R));
  flux_Momy = flux_Momy - (C * 0.5 * (rho_L * vy_L - rho_R * vy_R));
  flux_Energy = flux_Energy - (C * 0.5 * (en_L - en_R));

  return {flux_Mass, flux_Momx, flux_Momy, flux_Energy};
}

std::vector<double> linspace(double start, double stop, size_t num) {
  std::vector<double> result;

  double delta = (stop - start) / (num - 1);

  double curr = start;

  double count = 0;

  while (count < num) {
    result.push_back(curr);
    curr += delta;
    count++;
  }

  return result;
}

std::array<mtx, 2> meshgrid(const std::vector<double> &xv,
                            const std::vector<double> &yv) {
  mtx x_res(yv.size(), xv.size());
  mtx y_res(yv.size(), xv.size());

  for (size_t y = 0; y < x_res.rows(); ++y) {
    for (size_t x = 0; x < x_res.cols(); ++x) {
      x_res(y, x) = xv[x];
    }
  }

  for (size_t x = 0; x < y_res.cols(); ++x) {
    for (size_t y = 0; y < y_res.rows(); ++y) {
      y_res(y, x) = yv[y];
    }
  }

  return {x_res, y_res};
}

void simloop() {
  FILE *stream = fopen("output.bin", "wb");
  // TODO: Write N as first line of output.bin
  size_t N = 16;
  double boxsize = 1.0;
  double gamma = 5.0 / 3.0;
  double courant_fac = 0.4;
  double t = 0.0;
  double tEnd = 2.0;
  double tOut = 0.02;
  bool useSlopeLimiting = false;
  double plotRealTime = true;

  double dx = boxsize / (double)N;
  double vol = dx * dx;
  std::vector<double> xlin = linspace(0.5 * dx, boxsize - 0.5 * dx, N);
  auto [Y, X] = meshgrid(xlin, xlin);

  double w0 = 0.1;
  double sigma = 0.05 / std::sqrt(2.0);
  mtx rho = 1.0 + Matrix::filter_lt(Matrix::abs(Y - 0.5), 0.25);
  mtx vx = -0.5 + Matrix::filter_lt(Matrix::abs(Y - 0.5), 0.25);
  mtx vy = w0 * Matrix::sin(4 * M_PI * X) *
           (Matrix::exp(-1 * ((Y - 0.25) * (Y - 0.25)) / (2 * sigma * sigma)) +
            Matrix::exp(-1 * ((Y - 0.75) * (Y - 0.75)) / (2 * sigma * sigma)));
  mtx P = 2.5 * Matrix::ones(X.rows(), X.cols());

  auto [Mass, Momx, Momy, Energy] = getConserved(rho, vx, vy, P, gamma, vol);

  size_t outputCount = 1;

  while (t < tEnd) {
    auto [rho, vx, vy, P] = getPrimitive(Mass, Momx, Momy, Energy, gamma, vol);

    double dt =
        courant_fac * Matrix::min(dx / (Matrix::sqrt(gamma * P / rho) +
                                        Matrix::sqrt((vx * vx) + (vy * vy))));

    bool plotThisTurn = false;

    if (t + dt > outputCount * tOut) {
      dt = outputCount * tOut - t;
      plotThisTurn = true;
    }

    auto rho_d = getGradient(rho, dx);
    auto vx_d = getGradient(vx, dx);
    auto vy_d = getGradient(vy, dx);
    auto P_d = getGradient(P, dx);

    if (useSlopeLimiting) {
      rho_d = slopeLimit(rho, dx, rho_d[0], rho_d[1]);
      vx_d = slopeLimit(vx, dx, vx_d[0], vx_d[1]);
      vy_d = slopeLimit(vy, dx, vy_d[0], vy_d[1]);
      P_d = slopeLimit(P, dx, P_d[0], P_d[1]);
    }

    mtx rho_dx = std::move(rho_d[0]);
    mtx rho_dy = std::move(rho_d[1]);

    mtx vx_dx = std::move(vx_d[0]);
    mtx vx_dy = std::move(vx_d[1]);

    mtx vy_dx = std::move(vy_d[0]);
    mtx vy_dy = std::move(vy_d[1]);

    mtx P_dx = std::move(P_d[0]);
    mtx P_dy = std::move(P_d[1]);

    mtx rho_prime =
        rho -
        0.5 * dt * (vx * rho_dx + rho * vx_dx + vy * rho_dy + rho * vy_dy);

    mtx vx_prime = vx - 0.5 * dt * (vx * vx_dx + vy * vx_dy + (1 / rho) * P_dx);

    mtx vy_prime = vy - 0.5 * dt * (vx * vy_dx + vy * vy_dy + (1 / rho) * P_dy);

    mtx P_prime =
        P - 0.5 * dt * (gamma * P * (vx_dx + vy_dy) + vx * P_dx + vy * P_dy);

    auto [rho_XL, rho_XR, rho_YL, rho_YR] =
        extrapolateInSpaceToFace(rho_prime, rho_dx, rho_dy, dx);

    auto [vx_XL, vx_XR, vx_YL, vx_YR] =
        extrapolateInSpaceToFace(vx_prime, vx_dx, vx_dy, dx);
    auto [vy_XL, vy_XR, vy_YL, vy_YR] =
        extrapolateInSpaceToFace(vy_prime, vy_dx, vy_dy, dx);
    auto [P_XL, P_XR, P_YL, P_YR] =
        extrapolateInSpaceToFace(P_prime, P_dx, P_dy, dx);

    auto [flux_Mass_X, flux_Momx_X, flux_Momy_X, flux_Energy_X] =
        getFlux(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, gamma);
    auto [flux_Mass_Y, flux_Momy_Y, flux_Momx_Y, flux_Energy_Y] =
        getFlux(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, gamma);

    Mass = applyFluxes(Mass, flux_Mass_X, flux_Mass_Y, dx, dt);
    Momx = applyFluxes(Momx, flux_Momx_X, flux_Momx_Y, dx, dt);
    Momy = applyFluxes(Momy, flux_Momy_X, flux_Momy_Y, dx, dt);
    Energy = applyFluxes(Energy, flux_Energy_X, flux_Energy_Y, dx, dt);

    t += dt;

    if (plotThisTurn || t >= tEnd) {
      outputCount += 1;
      rho.writeToFile(stream);
    }
  }

  fclose(stream);
  printf("done!, outputcount: %zu\n", outputCount);
}
