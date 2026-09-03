// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#ifndef PISM_DEBRIS_TRANSPORT_HH
#define PISM_DEBRIS_TRANSPORT_HH

#include <memory>

#include "pism/coupler/DebrisModel.hh"
#include "pism/coupler/debris/DebrisInput.hh"
#include "pism/coupler/debris/EnglacialTransport.hh"
#include "pism/coupler/debris/GravitationalTransport.hh"
#include "pism/coupler/debris/SupraglacialTransport.hh"
#include "pism/coupler/debris/TerminusRemoval.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace debris {

//! @brief Prognostic englacial and supraglacial debris transport (Verhaegen and
//! Huybrechts, 2026).
/*!
  Debris arriving from the surrounding terrain (`DebrisInput`) is buried in the
  accumulation zone and advected englacially (`EnglacialTransport`); in the ablation zone
  it melts out and joins the supraglacial layer, which is advected by the surface velocity
  (`SupraglacialTransport`), redistributed downslope (`GravitationalTransport`) and
  removed at the glacier margin (`TerminusRemoval`).

  State: `debris_thickness` (m) and `englacial_debris_concentration` (kg m^-3). The model
  keeps a mass budget (inputs, melt-out, off-glacier output, losses) reported as scalar
  diagnostics; it does not rescale the debris mass.

  Select with `-debris transport`.
*/
class DebrisTransport : public DebrisModel {
public:
  DebrisTransport(std::shared_ptr<const Grid> grid);
  virtual ~DebrisTransport() = default;

  // sub-models and per-step fields (for diagnostics and tests)
  const EnglacialTransport &englacial() const;
  const SupraglacialTransport &supraglacial() const;
  const GravitationalTransport &gravitational() const;
  const TerminusRemoval &terminus() const;
  const DebrisInput &input() const;

  //! Englacial debris concentration (kg m^-3) at the end of the last step.
  const array::Array3D &concentration() const;

  //! Cumulative debris input (m of solid debris) since the start of the run.
  const array::Scalar &cumulative_input() const;
  //! Cumulative melt-out (m of solid debris) since the start of the run.
  const array::Scalar &cumulative_melt_out() const;
  //! Cumulative removal into the foreland (m) since the start of the run.
  const array::Scalar &cumulative_removal() const;

  //! `(1 - phi) rho`.
  double solid_density() const;
  //! Coefficient `C` of the debris-covered area fraction `1 - exp(-C h_d)`.
  double cover_fraction_coefficient() const;

  // mass budget (kg; reduced across processors)
  double englacial_mass() const;
  double supraglacial_mass() const;
  double total_mass() const;
  double initial_mass() const;
  double input_last_step() const;
  double melt_out_last_step() const;
  double output_last_step() const;
  double lost_last_step() const;
  double cumulative_input_mass() const;
  double cumulative_output_mass() const;
  double cumulative_lost_mass() const;
  //! `(M_englacial + M_supraglacial) - M_initial - input + output + lost` (kg)
  double conservation_error() const;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Inputs &inputs, double t, double dt);

  MaxTimestep max_timestep_impl(double t, const CFLData *cfl_data) const;

  const array::Scalar &debris_impl() const;

  std::set<VariableMetadata> state_impl() const;
  void write_state_impl(const OutputFile &output) const;

  DiagnosticList spatial_diagnostics_impl() const;
  TSDiagnosticList scalar_diagnostics_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;
  void init_timeseries_impl(const std::vector<double> &ts) const;
  void debris_time_series_impl(int i, int j, std::vector<double> &result) const;

private:
  void update_budget(const Geometry &geometry, double dt);

  double m_solid_density, m_vertical_cfl_ratio, m_cover_coefficient;
  bool m_do_gravity, m_do_terminus;

  // state
  array::Scalar2 m_debris_thickness;
  array::Array3D m_concentration;

  DebrisInput m_input;
  EnglacialTransport m_englacial;
  SupraglacialTransport m_supraglacial;
  GravitationalTransport m_gravity;
  TerminusRemoval m_terminus;

  // geometry the velocity field corresponds to (end of the previous step)
  array::Scalar1 m_H_old;
  array::CellType1 m_cell_type_old;

  // per-step contributions to the supraglacial layer (m) and cumulative amounts (m)
  array::Scalar m_surface_input, m_surface_loss;
  array::Scalar m_cumulative_input, m_cumulative_melt_out, m_cumulative_removal;

  // mass budget (kg)
  double m_mass_englacial, m_mass_supraglacial, m_mass_initial;
  double m_step_input, m_step_melt_out, m_step_output, m_step_lost;
  double m_cumulative_input_mass, m_cumulative_output_mass, m_cumulative_lost_mass;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_TRANSPORT_HH */
