#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdio.h>

/**
 * Calculates the conserved variables from the primitive variables
 *
 * @param rho - Density
 * @param vx - Velocity in the x direction
 * @param vy - Velocity in the y direction
 * @param P - Pressure
 * @param gamma - Ideal gas gamma
 * @param vol - Volume of the cell
 * @param Mass_o - Pointer to the buffer to store the mass
 * @param Momx_o - Pointer to the buffer to store the x momentum
 * @param Momy_o - Pointer to the buffer to store the y momentum
 * @param Energy_o - Pointer to the buffer to store the ernergy
 */
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

/**
 * Calculates the primitive variables from the conserved variables
 *
 * @param Mass - Mass
 * @param Momx - Momentum in the x direction
 * @param Momy - Momentum in the y direction
 * @param Energy - Energy
 * @param gamma - Ideal gas gamma
 * @param vol - Volume of the cell
 * @param rho_o - Pointer to the buffer to store the density
 * @param vx_o - Pointer to the buffer to store the x velocity
 * @param vy_o - Pointer to the buffer to store the y velocity
 * @param P_o - Pointer to the buffer to store the pressure
*/
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

/**
 * Calculate the gradients of an element in a field
 *
 * @param f - A variable in the field
 * @param f_down - The variable below f in the field
 * @param f_up - The variable above f in the field
 * @param f_right - The variable to the right of f in the field
 * @param f_left - The variable to the left of f in the field
 * @param dx - Size of a cell in the field
 * @param f_dx_o - Where to store the derivative of f in the x-direction
 * @param f_dy_o - Where to store the derivative of f in the y-direction 
 */
void getGradient(double f, double f_down, double f_up, double f_right, double f_left, double dx, double *f_dx_o, double *f_dy_o) {
  double f_dx = (f_up - f_down) / (2 * dx);
  double f_dy = (f_left - f_right) / (2 * dx);

  *f_dx_o = f_dx;
  *f_dy_o = f_dy;
}

/**
 * Performed spatial extrapolation on a cell in a field to each of 
 * the 4 faces of a cell (Its four neighbors up, down, left and right)
 *
 * @param f - Element in a field
 * @param f_dx - Derivative of f along the x-axis
 * @param f_dy - Derivate of f along the y-axis
 * @param dx - Cell size
 * @param f_XL_o - Pointer to the buffer to store the extrapolated variable to the left of f
 * @param f_XR_o - Pointer to the buffer to store the extrapolated variable to the right of f
 * @param f_YL_o - Pointer to the buffer to store the extrapolated variable below f
 * @param f_YR_o - Pointer to the buffer to store the extrapolated variable above f
*/
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

/**
 * Applies the fluxes to the conserved variables
 *
 * @param F - Conserved variable
 * @param flux_F_x - Flux of the conserved variable in the x direction
 * @param flux_F_x_up - Flux of the conserved variable above F in the x direction
 * @param flux_F_y_left - Flux of the conserved variable to the left of F in the y direction
 * @param flux_F_y - Flux of the conserved variable in the y direction
 * @param dx - Cell size
 * @param dt - Time step
 */
double applyFluxes(double F, double flux_F_x, double flux_F_x_up, double flux_F_y_left, double flux_F_y, double dx, double dt) {
  return F + (dt * dx) * (-flux_F_x + flux_F_x_up - flux_F_y + flux_F_y_left);
}

/**
 * Calculates the fluxes of a between two states
 *
 * @param rho_L - Density of the left state
 * @param rho_R - Density of the right state
 * @param vx_L - Velocity in the x direction of the left state
 * @param vx_R - Velocity in the x direction of the right state
 * @param vy_L - Velocity in the y direction of the left state
 * @param vy_R - Velocity in the y direction of the right state
 * @param P_L - Pressure of the left state
 * @param P_R - Pressure of the right state
 * @param gamma - Ideal gas gamma
 * @param flux_Mass_o - Pointer to the buffer to store the mass flux
 * @param flux_Momx_o - Pointer to the buffer to store the x momentum flux
 * @param flux_Momy_o - Pointer to the buffer to store the y momentum flux
 * @param flux_Energy_o - Pointer to the buffer to store the energy flux
*/
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

/**
 * Implementation of the numpy function linspace
 * generate a vector of num evenly distributed elements from 
 * start to stop
 *
 * @param start - Start of the vector
 * @param stop - End of the vector
 * @param num - Number of elements in the vector
 * @param output - Buffer to store the vector
 */
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

/**
 * Implementation of the numpy function meshgrid
 * generate a two square matricies by repeating the
 * row vector xv and the column vector yv
 *
 * @param xv - Row vector to generate x_res from
 * @param yv - Column vector to generate y_res from
 * @param N - Length of xv and yv vectors
 * @param x_res - Preallocated buffer to store matrix generated from xv in
 * @param y_res - Preallocated buffer to store matrix generated from yv in
*/
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

/**
 * Prints a matrix.
 *
 * @param a - buffer were the matrix is stored
 * @param w - Width of the matrix stored in a 
 * @param h - Height of the matrix stored in a 
 */
void printM(double *a, int w, int h) {
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      printf("%lf ", a[y * w + x]);
    }
    printf("\n");
  }
  printf("\n");
}

/**
 * Writes a matrix to a file in binary format.
 * The first 4 bytes in the file contain the matrixes width
 * The next 4 bytes contain its height
 * The following w * h sections of 8 bytes contain the values of 
 * the matrix
 *
 * @param stream - File output stream opened in "wb" mode to write to
 * @param m - buffer were matrix is stored 
 * @param w - Width of the matrix stored in m 
 * @param h - Height of the matrix stored in m 
 */
void writeToFile(FILE* stream, double * m, int w, int h) {
  fwrite(&w, sizeof(int), 1, stream); 
  fwrite(&h, sizeof(int), 1, stream); 
  fwrite(m, sizeof(double), w * h, stream);
}
