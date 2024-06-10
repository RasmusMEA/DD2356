
#ifndef FINITEVOL_H
#define FINITEVOL_H
#include <stdio.h> 

void getConserved(double rho, double vx, double vy, double P, double gamma, double vol, double *Mass_o, double *Momx_o, double *Momy_o,
                  double *Energy_o);

void getPrimitive(double Mass, double Momx, double Momy, double Energy, double gamma, double vol, double *rho_o, double *vx_o, double *vy_o,
                  double *P_o);

void getGradient(double f, double f_down, double f_up, double f_right, double f_left, double dx, double *f_dx_o, double *f_dy_o);

void extrapolateInSpaceToFace(double f, double f_dx, double f_dy, double dx, double *f_XL_o, double *f_XR_o, double *f_YL_o, double *f_YR_o);

double applyFluxes(double F, double flux_F_x, double flux_F_x_up, double flux_F_y_left, double flux_F_y, double dx, double dt);

void getFlux(double rho_L, double rho_R, double vx_L, double vx_R, double vy_L, double vy_R, double P_L, double P_R, double gamma,
             double *flux_Mass_o, double *flux_Momx_o, double *flux_Momy_o, double *flux_Energy_o); 

void linspace(double start, double stop, size_t num, double *output);

void meshgrid(double *xv, double *yv, size_t N, double *x_res, double *y_res);


void printM(double *a, int w, int h);
void writeToFile(FILE* stream, double * m, int w, int h);

#endif
