// Copyright (C) 2009-2010 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

#ifndef __enthalpyConverter_hh
#define __enthalpyConverter_hh

#include "NCVariable.hh"

//! Converts ice specific enthalpy to-and-from temperature and water content.
/*!
ONCE THIS CLASS IS CONFIRMED TO BE CORRECT, IT *CAN* GO IN src/base/pism_const.{hh|cc}

Use this way, for example within IceModel with NCConfigVariable config member:
\code
  #include "enthalpyConverter.hh"
  EnthalpyConverter EC(&config);  // this runs constructor, so make sure this
                                  //   is done after the creation/initialization of
                                  //   NCConfigVariable config
  ...
  for (...) {
    ...
    E_s = EC.getEnthalpyCTS(p);
    ... etc ...
  }   
\endcode

Some methods are functions which return the computed values or a boolean.  These
do no error checking.  Others return PetscErrorCode (and set arguments).  These
check either that the enthalpy is below that of liquid water, or that
the temperature (in K) is positive.

Specifically, getAbsTemp() gives return value 1 if the input enthalpy exceeded
that of liquid water, but puts the temperature of maximum-liquid-content 
temperate ice into its computed value for T.  And getWaterFraction() gives
return value of 1 under the same condition, but puts the maximum-liquid-content
into its computed value for omega.

The three methods that get the enthalpy from temperatures and liquid fractions, 
namely getEnth, getEnthPermissive, getEnthAtWaterFraction, are more strict about
error checking.  They call SETERRQ() if there arguments are invalid.
*/
class EnthalpyConverter {
public:
  EnthalpyConverter(const NCConfigVariable &config);

  PetscErrorCode viewConstants(PetscViewer viewer) const;

  double         getPressureFromDepth(double depth) const;
  double         getMeltingTemp(double p) const;
  double         getEnthalpyCTS(double p) const;
  PetscErrorCode getEnthalpyInterval(double p, double &E_s, double &E_l) const;
  double         getCTS(double E, double p) const;

  bool           isTemperate(double E, double p) const;
  bool           isLiquified(double E, double p) const;

  PetscErrorCode getAbsTemp(double E, double p, double &T) const;
  PetscErrorCode getPATemp(double E, double p, double &T_pa) const;

  PetscErrorCode getWaterFraction(double E, double p, double &omega) const;
  double         getWaterFractionLimited(double E, double p) const;

  PetscErrorCode getEnth(double T, double omega, double p, double &E) const;
  PetscErrorCode getEnthPermissive(double T, double omega, double p, double &E) const;
  PetscErrorCode getEnthAtWaterFraction(double omega, double p, double &E) const;

protected:
  double T_0, L, c_i, rho_i, g, p_air, beta;
};

#endif

