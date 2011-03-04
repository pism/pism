// Copyright (C) 2009-2011 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

//! Converts between specific enthalpy and temperature or liquid content.
/*!
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
namely getEnth(), getEnthPermissive(), getEnthAtWaterFraction(), are more strict about
error checking.  They call SETERRQ() if there arguments are invalid.
*/
class EnthalpyConverter {
public:
  EnthalpyConverter(const NCConfigVariable &config);
  virtual ~EnthalpyConverter() {}

  virtual PetscErrorCode viewConstants(PetscViewer viewer) const;

  virtual double         getPressureFromDepth(double depth) const;
  virtual double         getMeltingTemp(double p) const;
  virtual double         getEnthalpyCTS(double p) const;
  virtual PetscErrorCode getEnthalpyInterval(double p, double &E_s, double &E_l) const;
  virtual double         getCTS(double E, double p) const;

  virtual bool           isTemperate(double E, double p) const;
  virtual bool           isLiquified(double E, double p) const;

  virtual PetscErrorCode getAbsTemp(double E, double p, double &T) const;
  virtual PetscErrorCode getPATemp(double E, double p, double &T_pa) const;

  virtual PetscErrorCode getWaterFraction(double E, double p, double &omega) const;

  virtual PetscErrorCode getEnth(double T, double omega, double p, double &E) const;
  virtual PetscErrorCode getEnthPermissive(double T, double omega, double p, double &E) const;
  virtual PetscErrorCode getEnthAtWaterFraction(double omega, double p, double &E) const;

protected:
  double T_triple, L, c_i, rho_i, g, p_air, beta, T_tol;
  double T_0;
  bool   do_cold_ice_methods;
};


//! An EnthalpyConverter for use in verification tests.
/*!
Treats ice at any temperature as cold (= zero liquid fraction).  Makes absolute
temperature (in K) and enthalpy strictly proportional.  The reference
temperature is \f$T_0=0\f$:  \f$E = c_i T\f$.  The pressure-melting temperature,
and the enthalpy of the CTS (\f$E_s(p)\f$) is never reached because it is set
very, very high.

Note: Any instance of IceFlowLaw uses an EnthalpyConverter, and it is this
one when in verification mode.
 */
class ICMEnthalpyConverter : public EnthalpyConverter {
public:
  ICMEnthalpyConverter(const NCConfigVariable &config) : EnthalpyConverter(config) {
    T_0 = 0.0;
    // T_triple = 1.0e30;  // unreachable
    do_cold_ice_methods = true;
    // FIXME:  it *might* be nice to set these as overrides (?), but we have a "const"
    //   reference for config, and IceFlowLaw creates an EnthalpyConverter and
    //   *it* has a const reference ... so not for now
    //config.set("enthalpy_converter_reference_temperature",T_0);
    //config.set("water_triple_point_temperature",T_triple);
    //config.set_flag("do_cold_ice_methods",true);
  }

  virtual ~ICMEnthalpyConverter() {}

  /*! */
  virtual double getMeltingTemp(double /*p*/) const { return T_triple; }

  /*! */
  virtual PetscErrorCode getAbsTemp(double E, double /*p*/,
                                    double &T) const {
    T = E / c_i; return 0; }

  /*! */
  virtual PetscErrorCode getWaterFraction(double /*E*/, double /*p*/,
                                          double &omega) const {
    omega = 0.0; return 0; }

  /*! */
  virtual PetscErrorCode getEnth(double T, double /*omega*/, double /*p*/, 
                                 double &E) const {
    E = c_i * T; return 0; }

  /*! */
  virtual PetscErrorCode getEnthPermissive(double T, double /*omega*/, double /*p*/,
                                           double &E) const {
    E = c_i * T; return 0; }

  /*! */
  virtual PetscErrorCode getEnthAtWaterFraction(double /*omega*/, double p,
                                                double &E) const {
    E = getEnthalpyCTS(p); return 0; }

  virtual bool isTemperate(double /*E*/, double /*p*/) const { return false; }
};


#endif // __enthalpyConverter_hh

