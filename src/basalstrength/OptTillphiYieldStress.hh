// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#ifndef _PISMOPTTILLPHIYIELDSTRESS_H_
#define _PISMOPTTILLPHIYIELDSTRESS_H_

#include "MohrCoulombYieldStress.hh"

#include "pism/util/iceModelVec.hh"

namespace pism {

class IceModelVec2CellType;

//! @brief PISM's iteratively optimized basal yield stress model which applies the
//! Mohr-Coulomb model of deformable, pressurized till, with adjusted till friction angle.
class OptTillphiYieldStress : public MohrCoulombYieldStress {
public:
  OptTillphiYieldStress(IceGrid::ConstPtr g);
  virtual ~OptTillphiYieldStress() = default;

private:

  DiagnosticList diagnostics_impl() const;

  void iterative_phi_step(const IceModelVec2S &ice_surface_elevation,
                          const IceModelVec2S &bed_topography,
                          const IceModelVec2CellType &mask);

protected:

  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);

  void update_impl(const YieldStressInputs &inputs, double t, double dt);

  void restart_impl(const File &input_file, int record);
  //void init_impl(const YieldStressInputs &inputs);

  MaxTimestep max_timestep_impl(double t) const;

  IceModelVec2S m_diff_mask;
  IceModelVec2S m_diff_usurf;
  IceModelVec2S m_target_usurf;
  IceModelVec2S m_usurf;

  double m_last_time;
  double m_last_inverse_time;
  double m_dt_phi_inv;
};

} // end of namespace pism

#endif /* _PISMOPTTILLPHIYIELDSTRESS_H_ */
