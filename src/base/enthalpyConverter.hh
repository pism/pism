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

#ifndef __enthalpyConverter_hh
#define __enthalpyConverter_hh

#include "NCVariable.hh"

//! Converts enthalpy to-and-from temperature and water content.
/*!
THERE IS NO LONGER AN ATTEMPT TO INLINE THINGS.  IT WAS PREMATURE 
OPTIMIZATION ...

ONCE THIS CLASS IS CONFIRMED TO BE CORRECT, IT *CAN* GO IN src/base/pism_const.{hh|cc}

Use this way, for example within IceModel with NCConfigVariable config member:
\code
  #include "enthalpyConverter.hh"
  EnthalpyConverter EC(&config);  // this runs constructor, so make sure this
                                  //   is done after the creation/initialization of
                                  //   IceModel::config
  ...
  for (...) {
    ...
    E_s = EC.getEnthalpyCTS(p);
    ... etc ...
  }   
\endcode
*/
class EnthalpyConverter {
public:
  EnthalpyConverter(NCConfigVariable *config);

  PetscErrorCode viewConstants(PetscViewer viewer) const;

  // ice conversion methods
  double getPressureFromDepth(double depth) const;
  double getMeltingTemp(double p) const;
  double getEnthalpyCTS(double p) const;
  void   getEnthalpyInterval(double p, double &E_s, double &E_l) const;
  double getAbsTemp(double E, double p) const;
  double getPATemp(double E, double p) const;
  double getWaterFraction(double E, double p) const;
  double getEnth(double T, double omega, double p) const;
  double getEnthPermissive(double T, double omega, double p) const;
  double getCTS(double E, double p) const;

  // bedrock conversion methods
  double getEnthBedrock(double T) const;
  double getAbsTempBedrock(double E) const;

protected:
  double T_0, L, c_i, c_b, rho_i, g, p_air, beta;
};

#endif

