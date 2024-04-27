#include <array>
#include <tuple>
#include <vector>

#include "matrix.h"

#define mtx Matrix

using std::tuple;

std::array<mtx, 4> getConserved(const mtx &rho, const mtx &vx, const mtx &vy,
                                const mtx &P, double gamma, double vol);

std::array<mtx, 4> getPrimitive(const mtx &Mass, const mtx &Momx,
                                const mtx &Momy, const mtx &Energy,
                                double gamma, double vol);

std::array<mtx, 2> getGradient(const mtx &f, double dx);

std::array<mtx, 2> slopeLimit(const mtx &f, double dx, mtx f_dx, mtx f_dy);

std::array<mtx, 4> extrapolateInSpaceToFace(const mtx &f, const mtx &f_dx,
                                            const mtx &f_dy, double dx);

mtx applyFluxes(mtx F, const mtx &flux_F_x, const mtx &flux_F_y, double dx,
                double dt);

std::array<mtx, 4> getFlux(const mtx &rho_L, const mtx &rho_R, const mtx &vx_L,
                           const mtx &vx_R, const mtx &vy_L, const mtx &vy_R,
                           const mtx &P_L, const mtx &P_R, double gamma);

void simloop();

std::vector<double> linspace(double start, double stop, size_t num);

std::array<mtx, 2> meshgrid(const std::vector<double> &x,
                            const std::vector<double> &y);
