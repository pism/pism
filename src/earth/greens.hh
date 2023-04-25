// Copyright (C) 2007--2009, 2014, 2015, 2017, 2019 Ed Bueler
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

#ifndef __greens_hh
#define __greens_hh

#include <gsl/gsl_interp.h>  // for gsl_interp_accel
#include <gsl/gsl_spline.h>  // for gsl_spline

namespace pism {
namespace bed {

//! @brief The integrand in defining the elastic Green's function from
//! the [@ref Farrell] earth model.
/*!
 * For G^E(r), the Green's function of spherical layered elastic earth
 * model. From data in \ref LingleClark. See also \ref BLKfastearth.
 */
double ge_integrand(unsigned ndim, const double* args, void* data);

class greens_elastic {
public:
  greens_elastic();
  ~greens_elastic();
  double operator()(double r);
private:
  static const int N = 42;
  static const double rmkm[N];
  static const double GE[N];

  gsl_interp_accel* acc;
  gsl_spline*       spline;

  // disable copy constructor and the assignment operator:
  greens_elastic(const greens_elastic &other);
  greens_elastic& operator=(const greens_elastic&);
};

//! @brief Parameters used to access elastic Green's function from the
//! [@ref Farrell] earth model.
struct ge_data {
  double dx, dy;
  int    p, q;
  greens_elastic *G;
};

//! @brief Actually compute the response of the viscous half-space
//! model in \ref LingleClark, to a disc load.
double viscDisc(double t, double H0, double R0, double r,
                double rho, double rho_ice, double grav, double D, double eta);

} // end of namespace bed
} // end of namespace pism

#endif  /* __greens_hh */
