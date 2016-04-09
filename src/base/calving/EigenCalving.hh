/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef _PISMEIGENCALVING_H_
#define _PISMEIGENCALVING_H_

#include "base/util/iceModelVec.hh"
#include "base/util/PISMComponent.hh"

namespace pism {
namespace stressbalance {
class StressBalance;
}

class IceModelVec2CellType;

namespace calving {

class EigenCalving : public Component
{
public:
  EigenCalving(IceGrid::ConstPtr g, stressbalance::StressBalance *stress_balance);
  virtual ~EigenCalving();

  virtual void init();
  void update(double dt,
              double sea_level,
              const IceModelVec2S &bed_topography,
              IceModelVec2CellType &mask,
              IceModelVec2S &Href,
              IceModelVec2S &ice_thickness);

  MaxTimestep max_timestep();

  // empty methods that we're required to implement:
protected:
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO& nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                     IO_Type nctype);
  void compute_calving_rate(IceModelVec2CellType &mask);
  void update_strain_rates();

protected:
  IceModelVec2 m_strain_rates;
  IceModelVec2S m_thk_loss, m_horizontal_calving_rate;
  const int m_stencil_width;
  stressbalance::StressBalance *m_stress_balance;
  double m_K;
  bool m_restrict_timestep;
};

void remove_narrow_tongues(const IceModelVec2CellType &mask, IceModelVec2S &ice_thickness);

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMEIGENCALVING_H_ */
