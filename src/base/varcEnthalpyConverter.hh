// Copyright (C) 2011 Ed Bueler
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

#ifndef __varcEnthalpyConverter_hh
#define __varcEnthalpyConverter_hh

#include "enthalpyConverter.hh"

//! Enthalpy converter based on specific heat which is linear in temperature.
/*!
See equation (4.39) in \ref GreveBlatter2009,
        \f[C_i(T) = 146.3 + 7.253 T = 2009.0 + 7.253 (T - T_r)\f]
where \f$T\f$ is in Kelvin and the reference temperature is \f$T_r = 256.82\f$ K. 
 */
class varcEnthalpyConverter : public EnthalpyConverter {
public:
  varcEnthalpyConverter(const NCConfigVariable &config) : EnthalpyConverter(config) {}
  virtual ~varcEnthalpyConverter() {}

  virtual double         getEnthalpyCTS(double p) const;

  virtual PetscErrorCode getAbsTemp(double E, double p, double &T) const;

  virtual PetscErrorCode getEnth(double T, double omega, double p, double &E) const;

protected:
  double EfromT(double T) const;
  double TfromE(double E) const;
};

#endif // __varcEnthalpyConverter_hh

