// Copyright (C) 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _PO_RUNOFF_SMB_H_
#define _PO_RUNOFF_SMB_H_

#include "coupler/util/PScalarForcing.hh"
#include "POModifier.hh"

namespace pism {
namespace ocean {

/** "SMB Correction"
 *
 * This class implements the temperature-offset based SMB correction
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
class Runoff_SMB : public PScalarForcing<OceanModel,OceanModifier>
{
public:
  Runoff_SMB(IceGrid::ConstPtr g, OceanModel* in);
  virtual ~Runoff_SMB();

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void init_impl();
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result) const;

  double m_temp_to_runoff_a,
    m_temp_to_runoff_b,
    m_runoff_to_ocean_melt_a,
    m_runoff_to_ocean_melt_b,
    m_runoff_to_ocean_melt_power_alpha,
    m_runoff_to_ocean_melt_power_beta,
    m_current_forcing_0;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _PO_RUNOFF_SMB_H_ */
