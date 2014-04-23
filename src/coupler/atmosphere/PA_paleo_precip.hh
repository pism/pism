// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#ifndef _PA_PALEO_PRECIP_H_
#define _PA_PALEO_PRECIP_H_

#include "PScalarForcing.hh"
#include "PAModifier.hh"

namespace pism {

/** "Paleo-precipitation correction"
 *
 * This class implements the temperature-offset based precipitation
 * correction for use in paleo simulations.
 *
 * @f[
 * P_{G}(x,y,t) = P_{G}(x,y,0)exp\left[\frac{0.169}{d}\left({\Delta}T(t)+ {\Delta}T_{SC}(t)\right)\right]
 * @f]

 * where @f$P_G(x,y,0)@f$ is the precipitation for the present
 * Greenland ice sheet, obtained from the atmosphere model used as an
 * input of this class. The time dependent precipitation is
 * @f$P_G(x,y,t)@f$, @f$d@f$ is the @f${\delta}_{18}@f$ conversion factor
 * of @f$2.4^{\circ}C/\frac{0}{00}@f$, and @f${\Delta}T_{SC}@f$ is a
 * correction term for the change of altitude of the central dome
 * during the Greenland ice sheet's evolution. The coefficient
 * " @f$ 0.169/d @f$ " corresponds to a 7.3% change of precipitation rate
 * for every @f$1^{\circ}C@f$ of temperature change (Huybrechts et al.
 * 2002).
 *
 * One possible scheme for @f${\Delta}T_{SC}(t)@f$ (used in PISM) is
 * to take it to be zero, which regards the height correction as
 * belonging to the set of uncertainties related to the conversion
 * between isotopic and temperature signals.
 */
class PA_paleo_precip : public PScalarForcing<PISMAtmosphereModel,PAModifier>
{
public:
  PA_paleo_precip(IceGrid &g, const PISMConfig &conf, PISMAtmosphereModel* in);
  virtual ~PA_paleo_precip();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode init_timeseries(double *ts, unsigned int N);

  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);

  virtual PetscErrorCode precip_time_series(int i, int j, double *values);

  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);

protected:
  NCSpatialVariable air_temp, precipitation;
  double m_precipexpfactor;
  std::vector<double> m_scaling_values;
private:
  PetscErrorCode allocate_PA_paleo_precip();
};

} // end of namespace pism

#endif /* _PA_PALEO_PRECIP_H_ */
