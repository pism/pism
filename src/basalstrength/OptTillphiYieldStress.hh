// Copyright (C) 2011--2021 PISM Authors
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

//! Iterative optimization of the till friction angle.
class OptTillphiYieldStress : public MohrCoulombYieldStress {
public:
  OptTillphiYieldStress(IceGrid::ConstPtr g);
  virtual ~OptTillphiYieldStress() = default;

private:
  DiagnosticList diagnostics_impl() const;

  void update_tillphi(const IceModelVec2S &ice_surface_elevation,
                      const IceModelVec2S &bed_topography,
                      const IceModelVec2CellType &mask);

  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);
  void restart_impl(const File &input_file, int record);
  void update_impl(const YieldStressInputs &inputs, double t, double dt);

  MaxTimestep max_timestep_impl(double t) const;

  IceModelVec2S m_mask;
  IceModelVec2S m_usurf_difference;
  IceModelVec2S m_usurf_target;

  double m_last_inverse_time;
  double m_dt_phi_inv;

  double m_h_inv;
  double m_dhdt_conv;
  double m_dphi_max;
  double m_dphi_min;
  double m_phi_min;
  double m_phi_minup;
  double m_phi_max;
  double m_topg_min;
  double m_topg_max;
  double m_slope;
};

} // end of namespace pism

#endif /* _PISMOPTTILLPHIYIELDSTRESS_H_ */
