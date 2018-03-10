// Copyright (C) 2004-2007, 2015, 2017 Jed Brown and Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cassert>
#include <vector>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "matlablike.hh"

#include "greens.hh"

namespace pism {
namespace bed {

double ge_integrand(unsigned ndim, const double* xi_and_eta, void* parameters) {

  assert(ndim == 2);
  (void) ndim;

  // data here is from Lingle & Clark (1985), but claims to give G^E(r) for Farrell (1972)
  const double rmkm[42] =
    {0.0,    0.011,    0.111,  1.112,  2.224,  3.336,  4.448,  6.672,  8.896,  11.12,
     17.79,  22.24,    27.80,  33.36,  44.48,  55.60,  66.72,  88.96,  111.2,  133.4,
     177.9,  222.4,    278.0,  333.6,  444.8,  556.0,  667.2,  778.4,  889.6,  1001.0,
     1112.0, 1334.0,   1779.0, 2224.0, 2780.0, 3336.0, 4448.0, 5560.0, 6672.0, 7784.0,
     8896.0, 10008.0};
  // rm = rmkm * 1e3 (remember to convert to meters); GE /(10^12 rm) is
  //    vertical displacement in meters
  // (GE(r=0) has been computed by linear extrapolation:  GE(0) := -33.6488)
  const double GE[42] =
    {-33.6488, -33.64, -33.56, -32.75, -31.86, -30.98, -30.12, -28.44, -26.87, -25.41,
     -21.80,   -20.02, -18.36, -17.18, -15.71, -14.91, -14.41, -13.69, -13.01, -12.31,
     -10.95,   -9.757, -8.519, -7.533, -6.131, -5.237, -4.660, -4.272, -3.999, -3.798,
     -3.640,   -3.392, -2.999, -2.619, -2.103, -1.530, -0.292,  0.848,  1.676,  2.083,
      2.057,    1.643};

  struct ge_params* p = (struct ge_params*) parameters;

  const double
    dx        = p->dx,
    dy        = p->dy,
    xi_shift  = p->p * dx - xi_and_eta[0],
    eta_shift = p->q * dy - xi_and_eta[1],
    r         = sqrt(xi_shift * xi_shift + eta_shift * eta_shift);
  
  double z = 0.0;
  if (r < 0.01) {
    z = GE[0]/(rmkm[1] * 1.0e3 * 1.0e12);
  } else if (r > rmkm[41] * 1.0e3) {
    z = 0.0;
  } else {
    z = interp1_linear((double*)rmkm, (double*)GE, 42, r / 1.0e3) / (r * 1.0e12);
  }
  return z;
}


double vd_integrand (double kappa, void* parameters) {
  // Matlab:  function y=integrand(kappa,rg,D,t,eta,R0,rk)
  //            beta=rg + D*kappa.^4;
  //            expdiff=exp(-beta*t./(2*eta*kappa))-ones(size(kappa));
  //            y=expdiff.*besselj(1.0,kappa*R0).*besselj(0.0,kappa*rk)./beta;

  struct vd_params*  p = (struct vd_params*) parameters;

  const double
    t       = p->t,
    R0      = p->R0,
    rk      = p->rk,
    rho     = p->rho,
    grav    = p->grav,
    D       = p->D,
    eta     = p->eta,
    beta    = rho * grav + D * pow(kappa, 4.0),
    expdiff = exp(-beta * t / (2.0 * eta * kappa)) - 1.0;

  return expdiff * gsl_sf_bessel_J1(kappa * R0) * gsl_sf_bessel_J0(kappa * rk) / beta;
}

/*!
 * Compute the response of the viscous half-space model in [@ref LingleClark] to a disc load.
 *
 * @param[in] t time in seconds
 * @param[in] H0 thickness of the disc load, meters
 * @param[in] R0 radius of the disc load, meters
 * @param[in] r radius, meters
 * @param[in] rho mantle density, kg/m3
 * @param[in] rho_ice ice (load) density, kg/m3
 * @param[in] grav acceleration due to gravity, m/s2
 * @param[in] D lithosphere flexural rigidity, N meter
 * @param[in] eta mantle viscosity, Pascal second
 */
double viscDisc(double t, double H0, double R0, double r, 
                double rho, double rho_ice, double grav, double D, double eta) {
  // t in seconds; H0, R0, r in meters

  const double ABSTOL          = 1.0e-10;
  const double RELTOL          = 1.0e-14;
  const int    N_gsl_workspace = 1000;
  const int    lengthpts       = 142;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc(N_gsl_workspace);
  
  // Matlab:  pts=[10.^(-3:-0.05:-10) 1.0e-14];
  std::vector<double> pts(lengthpts);
  for (int j=0; j < lengthpts-1; j++) {
    pts[j] = pow(10.0, -3.0 - 0.05 * j);
  }
  pts[lengthpts-1] = 1.0e-14;

  // result=quadl(@integrand,pts(1),100.0*pts(1),TOL,0,rg,D,t,eta,R0,rk); % kap->infty tail
  gsl_function      F;
  struct vd_params  params = { t, R0, r, rho, grav, D, eta };
  double            result, error;
  F.function = &vd_integrand;
  F.params = &params;
  // regarding tolerance: request is for convergence of all digits and relative tolerance RELTOL
  gsl_integration_qag (&F, pts[1], 100.0*pts[1], ABSTOL, RELTOL, N_gsl_workspace, 
                       GSL_INTEG_GAUSS21, w, &result, &error);

  double sum = result;
  // for j=1:length(pts)-1
  //   result=result+quadl(@integrand,pts(j+1),pts(j),TOL,0,rg,D,t,eta,R0,rk);
  // end
  for (int j = 0; j < lengthpts - 1; j++) {
    gsl_integration_qag (&F, pts[j+1], pts[j], ABSTOL, RELTOL, N_gsl_workspace, 
                         GSL_INTEG_GAUSS21, w, &result, &error);
    sum += result;
  }
  
  gsl_integration_workspace_free(w);
  // u(k)=rhoi*g*H0*R0*result;
  return rho_ice * grav * H0 * R0 * sum;
}

} // end of namespace bed
} // end of namespace pism
