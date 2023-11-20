/*
   Copyright (C) 2011, 2013, 2014, 2016, 2023 Ed Bueler

   This file is part of PISM.

   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
   version.

   PISM is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with PISM; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "pism/verification/tests/exactTestO.h"

static const double beta_CC   = 7.9e-8; /* K Pa-1; Clausius-Clapeyron constant [@ref Luethi2002] */
static const double T_triple  = 273.15; /* K; triple point of pure water */
static const double L         = 3.34e5; /* J kg-1; latent heat of fusion for water [@ref AschwandenBlatter] */
static const double grav      = 9.81;   /* m/s^2; accel of gravity */
static const double rho_ICE   = 910.0;  /* kg/(m^3)  density of ice */
static const double k_ICE     = 2.10;   /* J/(m K s) = W/(m K)  thermal conductivity of ice */
static const double k_BEDROCK = 3.0;    /* J/(m K s) = W/(m K)  thermal conductivity of bedrock */
static const double H0        = 3000.0; /* m */
static const double Ts        = 223.15; /* K */
static const double G         = 0.042;  /* W/(m^2) = J m-2 s-1 */

/*! \brief Implements an exact solution for basal melt rate.  Utterly straightforward arithmetic. */
/*!
Assumes a steady-state temperature profile in ice and bedrock.  This steady
profile is driven by geothermal flux `G` from the crust below the bedrock
thermal layer.  The heat flux is everywhere upward in the bedrock thermal layer,
and it is equal by construction to `G`.  The heat flux upward in the ice is
also a constant everywhere in the ice, but its value is smaller than `G`.  This
imbalance is balanced by the basal melt rate `bmelt`.

Geometry and dynamics context:  The top of the ice is flat so the ice
does not flow; the ice thickness has constant value `H0`.  The ice surface has
temperature `Ts` which is below the melting point.  The basal melt rate
`bmelt` does not contribute to the vertical velocity in the ice.  The ice
pressure is hydrostatic: \f$p = \rho_i g (h-z)\f$.

The basic equation relates the fluxes in the basal ice, and in the top of the
bedrock layer, with the basal melt rate `bmelt` \f$= - M_b / \rho_i \f$.  The
equation is from [\ref AschwandenBuelerKhroulevBlatter],
  \f[ M_b H + (\mathbf{q} - \mathbf{q_{lith}}) \cdot \mathbf{n} = F_b + Q_b. \f]
Here \f$-M_b\f$ is the basal melt rate in units of mass per area per time.
In the above equation the basal friction heating
\f$F_b\f$ is zero and the subglacial aquifer enthalpy flux \f$Q_b\f$ includes no
horizontal velocity.  (Note that \f$Q_b\f$ is the heat delivered by subglacial
water to the base of the ice.)  We assume the subglacial water is at the ice
overburden pressure \f$p_0=\rho_i g H_0\f$, and we assume that the temperate
layer at the base of the ice has zero thickness, so \f$\omega = 0\f$.  Thus
  \f[ H_l(p_b) = H_l(p_0) = H_s(p_0) + L, \f]
  \f[ H = H_s(p_0) + \omega L = H_s(p_0), \f]
  \f[ Q_b = H_l(p_0) M_b. \f]

The basic equation therefore reduces to
  \f[ \mathbf{q} \cdot \mathbf{n} - G = M_b L \f]
or
  \f[ \text{\texttt{bmelt}} = -\frac{M_b}{\rho_i}
                = \frac{G - \mathbf{q} \cdot \mathbf{n}}{L \rho_i}. \f]
On the other hand, the temperature in the ice is the steady-state result wherein
the upward flux is constant and the (absolute) temperature is a linear function
of the elevation \f$z\f$, so
  \f[ \mathbf{q} \cdot \mathbf{n} = - k_i \frac{T_s - T_m(p)}{H_0}.\f]

The temperature in the ice (\f$0 \le z \le H_0\f$) is this linear function,
  \f[ T(z) = T_m(p) + \frac{T_s - T_m(p)}{H_0} z \f]
and in the bedrock (\f$z \le 0\f$) is also linear,
  \f[ T_b(z) = T_m(p) - \frac{G}{k_b} z. \f]

This method implements these formulas.  It should be called both when setting-up
a verification test by setting temperature at different elevations within
the ice and bedrock, and when doing the verification itself by checking against
the exact `bmelt` value.
 */
int exactO_old(double z, double *TT, double *Tm, double *qice, double *qbed, double *bmelt) {

  double P_base;

  P_base = rho_ICE * grav * H0;     /* Pa; hydrostatic approximation to pressure
                                       at base */

  *Tm = T_triple - beta_CC * P_base;/* K; pressure-melting point at base */

  *qice = - k_ICE * (Ts - *Tm) / H0;/* J m-1 K-1 s-1 K m-1 = J m-2 s-1 = W m-2;
                                       equilibrium heat flux everywhere in ice */

  *qbed = G;                        /* J m-2 s-1 = W m-2; equilibrium heat flux
                                       everywhere in bedrock */

  *bmelt = (G - *qice) / (L * rho_ICE);/* J m-2 s-1 / (J kg-1 kg m-3) = m s-1;
                                          ice-equivalent basal melt rate */

  if (z > H0) {
    *TT = Ts;                       /* K; in air above ice */
  } else if (z >= 0.0) {
    *TT = *Tm + (Ts - *Tm) * (z / H0); /* in ice */
  } else {
    *TT = *Tm - (G / k_BEDROCK) * z;   /* in bedrock */
  }
  return 0;
}

struct TestOParameters exactO(double z) {
  struct TestOParameters result;

  exactO_old(z, &result.TT, &result.Tm, &result.qice, &result.qbed, &result.bmelt);

  return result;
}
