#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

void linSpace(double from, double to, int number, double *output); // impl

void meshGrid(double *arr, int N, double *output1, double *output2); // impl

void init_rho(double *Y, int N, double *rho_o);

void init_vx(double *Y, int N, double *vx_o);

void init_vy(double *X, double *Y, double w0, double sigma, int N,
             double *vy_o);

void init_P(int N, double *P_o);

void getConserved(double *rho, double *vx, double *vy, double *P, int N,
                  double gamma, double vol, double *Mass_o, double *Momx_o,
                  double *Momy_o, double *Energy_o);

void getPrimitive(double *Mass, double *Momx, double *Momy, double *Energy,
                  int N, double gamma, double vol, double *rho_o, double *vx_o,
                  double *vy_o, double *P_o);

void getGradient(double *f, int N, double dx, double *f_dx_o, double *f_dy_o);

void getRhoPrime(double *rho, double *rho_dx, double *rho_dy, double *vx,
                 double *vx_dx, double *vy, double *vy_dy, double dt, int N,
                 double *rho_prime_o);

void getVxPrime(double *rho, double *vx, double *vx_dx, double *vx_dy,
                double *vy, double *P_dx, double dt, int N, double *vx_prime_o);

void getVyPrime(double *rho, double *vx, double *vy, double *vy_dx,
                double *vy_dy, double *P_dy, double dt, int N,
                double *vy_prime_o);

void getPPrime(double *vx, double *vx_dx, double *vy, double *vy_dy, double *P,
               double *P_dx, double *P_dy, double dt, double gamma, int N,
               double *P_prime_o);

void extrapolate(double *f, double *f_dx, double *f_dy, int N, double dx,
                 double *f_XL_o, double *f_XR_o, double *f_YL_o,
                 double *f_YR_o);

void getFlux(double *rho_L, double *rho_R, double *vx_L, double *vx_R,
             double *vy_L, double *vy_R, double *P_L, double *P_R, int N,
             double gamma, double *flux_Mass_o, double *flux_Momx_o,
             double *flux_Momy_o, double *flux_Energy_o);

void applyFluxes(double *F, double *flux_F_X, double *flux_F_Y, int N,
                 double dx, double dt, double *output);

double getDt(double *rho, double *vx, double *vy, double *P, int N,
             double courant_fac, double gamma, double dx);

#endif
