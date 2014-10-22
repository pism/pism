// Copyright (C) 2011, 2014 Ed Bueler
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

#ifndef __varcEnthalpyConverter_hh
#define __varcEnthalpyConverter_hh

#include "enthalpyConverter.hh"

namespace pism {

//! Enthalpy converter based on specific heat which is linear in temperature.
/*!
  This EnthalpyConverter object converts between enthalpy and temperature using
  linear-in-temperature specific heat capacity.  We implement only 
  Equation (4.39) in [\ref GreveBlatter2009],
  \f[ C(T) = 146.3 + 7.253 T = c_i + 7.253 (T - T_r) \f]
  where \f$T\f$ is in Kelvin, \f$c_i = 2009\,\, \text{J}\,\text{kg}^{-1}\,\text{K}^{-1}\f$,
  and the reference temperature is \f$T_r = 256.81786846822\f$ K.
*/
class varcEnthalpyConverter : public EnthalpyConverter {
public:
  varcEnthalpyConverter(const Config &config)
    : EnthalpyConverter(config),
      T_r(256.81786846822),
      c_gradient(7.253)
  {
  }
  virtual ~varcEnthalpyConverter() {}

  virtual double         getEnthalpyCTS(double p) const;

  virtual PetscErrorCode getAbsTemp(double E, double p, double &T) const;

  virtual PetscErrorCode getEnth(double T, double omega, double p, double &E) const;

  /*!
    Equation (4.39) in [\ref GreveBlatter2009] is
    \f$C(T) = c_i + 7.253 (T - T_r)\f$, with a reference temperature
    \f$T_r = 256.82\f$ K.
  */
  virtual double c_from_T(double T) const
  { return c_i + c_gradient * (T - T_r); }

  virtual double c_from_enth(double E, double p) const
  {
    double T;
    getAbsTemp(E, p, T);
    return c_from_T(T);
  }

protected:
  const double T_r,  //!< reference temperature in the parameterization of C(T)
    c_gradient;      //!< \brief the rate of change of C with respect to T in
  //!< the parameterization of C(T)
  double EfromT(double T) const;
  double TfromE(double E) const;
};

} // end of namespace pism

#endif // __varcEnthalpyConverter_hh

