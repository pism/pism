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
ONCE THESE ARE CONFIRMED TO BE CORRECT, THEY CAN GO IN src/base/pism_const.hh
*/

//! Get pressure in ice from depth below surface using the hydrostatic assumption.
/*! If \f$d\f$ is the depth, so generally \f$d = \f$H[i][j] - z[k], then
     \f[ p = p_{\text{air}}  + \rho_i * g * d. \f]
 */ 
inline PetscScalar getPressureFromDepth(NCConfigVariable &config, PetscScalar depth) {
  static const PetscScalar p_air = config.get("surface_pressure"); // Pa
  if (depth <= 0.0) { // at or above surface of ice
    return p_air;
  } else {
    static const PetscScalar g     = config.get("earth_gravity"),
                             rho_i = config.get("ice_density");
    return p_air + rho_i * g * depth;
  }
}


//! Get melting temperature from pressure p.
/*!
     \f[ T_m(p) = T_0 - \beta p. \f]
 */ 
inline PetscScalar getMeltingTemp(NCConfigVariable &config, PetscScalar p) {
  static const PetscScalar T_0  = config.get("water_melting_temperature"),    // K
                           beta = config.get("beta_CC");                      // K Pa-1
  return T_0 - beta * p;
}


//! Get enthalpy at phase transition endpoints, and pressure melting temperature, from pressure p.
/*! From \ref AschwandenBlatter2009,
     \f[ H_s(p) = -L + H_l(p), \f]
     \f[ H_l(p) = c_w T_m(p). \f]
 */ 
inline void getEnthalpyInterval(NCConfigVariable &config, PetscScalar p, 
                                PetscScalar &H_s, PetscScalar &H_l) {
  static const PetscScalar c_w  = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                           L    = config.get("water_latent_heat_fusion");     // J kg-1
  H_l = c_w * getMeltingTemp(config,p);
  H_s = - L + H_l;
}


//! Get absolute ice temperature (K) from enthalpy H and pressure p.
/*! From \ref AschwandenBlatter2009, equation (12)
     \f[ T=T(H,p) = \begin{cases} 
                       c_i^{-1} (H-H_s(p)) + T_m(p), & H < H_s(p), \\
                       T_m(p), &                       H_s(p) \le H < H_l(p).
                    \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$H \ge H_l(p)\f$.
 */
inline PetscScalar getAbsTemp(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  const PetscScalar T_m = getMeltingTemp(config,p);
  PetscScalar H_s, H_l;
  getEnthalpyInterval(config, p, H_s, H_l);
  if (H < H_s) {
    static const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
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
     \f[ T_{pa} = T_{pa}(H,p) = \begin{cases} 
                                  c_i^{-1} (H-H_s(p)) + T_0, & H < H_s(p), \\
                                  T_0,                       & H_s(p) \le H < H_l(p).
                                \end{cases} \f]
 */
inline PetscScalar getPATemp(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar H_s, H_l;
  getEnthalpyInterval(config, p, H_s, H_l);
  static const PetscScalar T_0  = config.get("water_melting_temperature");     // K
  if (H < H_s) {
    static const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
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
/*! From \ref AschwandenBlatter2009, equation (12),
   \f[ \omega=\omega(H,p) = \begin{cases}
                               0.0,            & H \le H_s(p), \\
                               (H-H_s(p)) / L, & H_s(p) < H < H_l(p).
                            \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$H \ge H_l(p)\f$.
 */
inline PetscScalar getWaterFraction(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar H_s, H_l;
  getEnthalpyInterval(config, p, H_s, H_l);
  if (H <= H_s) { // two cases in (12)
    return 0.0;
  } else if (H < H_l) {
    static const PetscScalar L = config.get("water_latent_heat_fusion");     // J kg-1
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
/*!  From \ref AschwandenBlatter2009,
   \f[ c = (1 - \omega) c_i + \omega c_w, \f]
   \f[ H = H_s(p) + c (T-T_0) + \omega L. \f]
Also does checks that \f$T\f$ is a reasonable ice temperature (i.e. it is at or below \f$T_0\f$),
and that \f$\omega\f$ is a meaningful fraction (i.e. it is in the interval \f$(0,1)\f$).

(Which is correct, ``\f$T-T_0\f$'' or ``\f$T-T_m\f$'' in this case.  Note that in this
case \f$T\f$ is the absolute, not the pressure-adjusted, temperature.)
 */
inline PetscScalar getEnth(NCConfigVariable &config, PetscScalar T, PetscScalar omega, PetscScalar p) {
  if ((omega < 0.0) || (1.0 < omega)) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getEnth(): water fraction omega not in range [0,1]; ending ... \n\n");
    PetscEnd();
  }
  static const PetscScalar T_0 = config.get("water_melting_temperature");    // K
  if (T > T_0 + 0.000001) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getEnth(): T exceeds T_0 so we have liquid water; ending ... \n\n");
    PetscEnd();
  }
  PetscScalar H_s, H_l;
  getEnthalpyInterval(config, p, H_s, H_l);
  static const PetscScalar c_w = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                           c_i = config.get("ice_specific_heat_capacity"),   // J kg-1 K-1
                           L = config.get("water_latent_heat_fusion");       // J kg-1
  static const PetscScalar c = (1.0 - omega) * c_i + omega * c_w;
  return H_s + c * (T - T_0) + omega * L;
}


//! Compute enthalpy from absolute temperature within the bedrock using conditions at z=0 (bottom of ice and top of the bedrock) to set enthalpy scale.
/*! If \f$c_b\f$ is the specific heat capacity of the bedrock, if \c Enth_zero = \f$H(z=0)\f$ is
the enthalpy at the bottom of the ice, and if \c Temp_zero = \f$T(z=0)\f$ is the temperature at
the top of the bedrock then we compute and return
\f[ H(z < 0) = H(z=0) + c_b (T(z<0) - T(z=0)) \f]
where \c Tb = \f$T(z<0)\f$.  Input temperatures are assumed to be absolute (not pressure-adjusted).
Generally \f$T(z<0) > T(z=0)\f$, in which case the resulting enthalpy
value comes out higher than \f$H(z=0)\f$, and it scales linearly with the increasing temperature
as we descend into the bedrock.
 */
inline PetscScalar getEnthBedrock(NCConfigVariable &config, 
                                  PetscScalar Enth_zero, PetscScalar Temp_zero, PetscScalar Tb) {
  static const PetscScalar c_b = config.get("bedrock_thermal_specific_heat_capacity");   // J kg-1 K-1
  return Enth_zero + c_b * (Tb - Temp_zero);
}


//! Inverse function from getEnthBedrock().
/*! In same notation as for getEnthBedrock(),
\f[ T(z < 0) = \frac{H(z<0) - H(z=0)}{c_b} + T(z=0) \f]
 */
inline PetscScalar getAbsTempBedrock(NCConfigVariable &config, 
                                     PetscScalar Enth_zero, PetscScalar Temp_zero, PetscScalar Enthb) {
  static const PetscScalar c_b = config.get("bedrock_thermal_specific_heat_capacity");   // J kg-1 K-1
  return ((Enthb - Enth_zero) / c_b) + Temp_zero;
}


#endif

