// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef __enthalpyHelper_hh
#define __enthalpyHelper_hh


/*
My intent about these is that, *once everything works*, they might get massaged a bit
to be more efficient.  For instance, getEnthalpyInterval() is being used in some places
just to compute the pressure-melting temperature T_m.  More broadly, the repeated use
of config.get() will have an efficiency cost, so these inline procedures could
have versions which take as arguments the things they currently get through config.get().
*/

static inline PetscScalar getPressureFromDepth(NCConfigVariable &config, PetscScalar depth) {
  const PetscScalar p_air = config.get("surface_pressure"); // Pa
  if (depth <= 0.0) { // at or above surface of ice
    return p_air;
  } else {
    const PetscScalar g     = config.get("earth_gravity"),
                      rho_i = config.get("ice_density");
    return p_air + rho_i * g * depth;
  }
}


//! Get enthalpy at phase transition endpoints, and pressure melting temperature, from pressure p.
/*! From \ref AschwandenBlatter2009,
     \f[ T_m(p) = T_0 - \beta p, \f]
     \f[ H_l(p) = c_w T_m(p), \f]
     \f[ H_s(p) = -L + H_l(p). \f]
Also returns H_s.
 */ 
static inline PetscScalar getEnthalpyInterval(NCConfigVariable &config, PetscScalar p, 
                                       PetscScalar &T_m, PetscScalar &H_l, PetscScalar &H_s) {
  const PetscScalar T_0  = config.get("water_melting_temperature"),    // K
                    beta = config.get("beta_CC"),                      // K Pa-1
                    c_w  = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                    L    = config.get("water_latent_heat_fusion");     // J kg-1
  T_m = T_0 - beta * p;
  H_l = c_w * T_m;
  H_s = - L + H_l;
  return H_s;
}


//! Get absolute ice temperature (K) from enthalpy H and pressure p.
/*! From \ref AschwandenBlatter2009,
     \f[ T=T(H,p) = \begin{cases} 
                       c_i^{-1} (H-H_s(p)) + T_m(p), & H < H_s(p), \\
                       T_m(p), &                       H_s(p) \le H < H_l(p).
                    \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$H \ge H_l(p)\f$.
  
See getEnthalpyInterval() for calculation of \f$T_m(p)\f$, \f$H_l(p)\f$, and \f$H_s(p)\f$. 
 */
static inline PetscScalar getAbsTemp(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  // implement T part of eqn (12) in AB2009, but bonk if liquid water
  if (H < H_s) {
    const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
    return ((H - H_s) / c_i) + T_m;
  } else if (H < H_l) { // two cases in (12)
    return T_m;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getAbsTemp():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return T_m;
  }
}


//! Get pressure-adjusted ice temperature (K) from enthalpy H and pressure p.
/*! See getAbsTemp(), and compare:
     \f[ T_pa = T_pa(H,p) = \begin{cases} 
                       c_i^{-1} (H-H_s(p)) + T_0, & H < H_s(p), \\
                       T_0,                       & H_s(p) \le H < H_l(p).
                    \end{cases} \f]
 */
static inline PetscScalar getPATemp(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  const PetscScalar T_0  = config.get("water_melting_temperature");     // K
  if (H < H_s) {
    const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
    return ((H - H_s) / c_i) + T_0;
  } else if (H < H_l) { // two cases in (12)
    return T_0;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getPATemp():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return T_0;
  }
}


//! Get liquid water fraction from enthalpy H and pressure p.
/*! From \ref AschwandenBlatter2009,
   \f[ \omega=\omega(H,p) = \begin{cases}
                               0.0,            & H \le H_s(p), \\
                               (H-H_s(p)) / L, & H_s(p) < H < H_l(p).
                            \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$H \ge H_l(p)\f$.
  
See getEnthalpyInterval() for calculation of \f$T_m(p)\f$, \f$H_l(p)\f$, and \f$H_s(p)\f$. 
 */
static inline PetscScalar getWaterFraction(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  // implement omega part of eqn (12) in AB2009, but bonk if liquid water
  if (H <= H_s) { // two cases in (12)
    return 0.0;
  } else if (H < H_l) {
    const PetscScalar L = config.get("water_latent_heat_fusion");     // J kg-1
    return (H - H_s) / L;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getWaterFraction():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return 1.0;
  }
}


//! Compute enthalpy from absolute temperature, liquid water fraction, and pressure p.
static inline PetscScalar getEnth(NCConfigVariable &config, PetscScalar T, PetscScalar omega, PetscScalar p) {
  if ((omega < 0.0) || (1.0 < omega)) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getEnth(): water fraction omega not in range [0,1]; ending ... \n\n");
    PetscEnd();
  }
  const PetscScalar T_0 = config.get("water_melting_temperature");    // K
  if (T > T_0 + 0.000001) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getEnth(): T exceeds T_0 so we have liquid water; ending ... \n\n");
    PetscEnd();
  }
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  const PetscScalar c_w = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                    c_i = config.get("ice_specific_heat_capacity"),   // J kg-1 K-1
                    L = config.get("water_latent_heat_fusion");       // J kg-1
  const PetscScalar c = (1.0 - omega) * c_i + omega * c_w;
  return H_s + c * (T - T_0) + omega * L;
}


#endif

