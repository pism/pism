// Copyright (C) 2011--2022 PISM Authors
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

namespace pism {

//! Iterative optimization of the till friction angle.
class OptTillphiYieldStress : public MohrCoulombYieldStress {
public:
  OptTillphiYieldStress(IceGrid::ConstPtr g);
  virtual ~OptTillphiYieldStress() = default;

private:
  DiagnosticList diagnostics_impl() const;

  void update_tillphi(const array::Scalar &ice_surface_elevation,
                      const array::Scalar &bed_topography,
                      const array::CellType &mask);

  void init_t_last(const File &input_file);
  void init_usurf_target(const File &input_file);

  void init_impl(const YieldStressInputs &inputs);
  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);
  void restart_impl(const File &input_file, int record);
  void update_impl(const YieldStressInputs &inputs, double t, double dt);

  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  array::Scalar1 m_mask;
  array::Scalar1 m_usurf_difference;
  array::Scalar1 m_usurf_target;

  double m_dphi_scale;

  // convergence threshold:
  double m_dhdt_min;

  // lower and upper bounds of a tillphi adjustment:
  double m_dphi_min;
  double m_dphi_max;

  // constants defining the lower bound of tillphi:
  double m_phi0_min;
  double m_phi0_max;
  double m_topg_min;
  double m_topg_max;

  // the upper bound of tillphi:
  double m_phi_max;

  //! time of the last till friction angle update
  double m_t_last;
  //! Update interval in seconds
  double m_update_interval;
  //! Temporal resolution to use when checking whether it's time to update
  double m_t_eps;
  //! Name of the variable used to store the last update time.
  std::string m_time_name;
};

} // end of namespace pism

#endif /* _PISMOPTTILLPHIYIELDSTRESS_H_ */
