// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 PISM Authors
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

#include "YieldStress.hh"
#include "MohrCoulombYieldStress.hh"

#include "pism/util/iceModelVec.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {

class IceModelVec2CellType;

//! @brief PISM's iteratively optimized basal yield stress model which applies the
//! Mohr-Coulomb model of deformable, pressurized till, with adjusted till friction angle.
class OptTillphiYieldStress : public MohrCoulombYieldStress {
public:
  OptTillphiYieldStress(IceGrid::ConstPtr g);
  virtual ~OptTillphiYieldStress();

  //void set_till_friction_angle(const IceModelVec2S &input);

private:

  DiagnosticList diagnostics_impl() const;

private:

  void iterative_phi_step(const IceModelVec2S &ice_surface_elevation,
                          const IceModelVec2S &bed_topography,
                          const IceModelVec2CellType &mask);

protected:

  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);

  void update_impl(const YieldStressInputs &inputs, double t, double dt);


  //void restart_impl(const File &input_file, int record);
  //void init_impl(const YieldStressInputs &inputs);

  //void define_model_state_impl(const File &output) const;
  //void write_model_state_impl(const File &output) const;

  //void finish_initialization(const YieldStressInputs &inputs);


  MaxTimestep max_timestep_impl(double t) const;

  bool m_iterative_phi;
  IceModelVec2S m_target_usurf,m_diff_usurf,m_usurf,m_diff_mask;
  double m_last_time,m_last_inverse_time,dt_phi_inv;

};

} // end of namespace pism

#endif /* _PISMOPTTILLPHIYIELDSTRESS_H_ */