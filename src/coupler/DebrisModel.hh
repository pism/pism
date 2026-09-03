// Copyright (C) 2026 Constantine Khroulev, Ricarda Winkelmann and Andy Aschwanden
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

#ifndef PISM_DEBRIS_MODEL_HH
#define PISM_DEBRIS_MODEL_HH

#include <vector>

#include "pism/util/Component.hh"

namespace pism {

//! @brief Supraglacial debris models and modifiers: provide the debris thickness
//! (and, through debris::IceMeltEnhancement, its effect on melt).
namespace debris {

//! Inputs of a debris model.
/*!
  Only `geometry` is required by every model; the prognostic transport model also needs
  the 3D velocity field and the surface mass balance applied during the current step.
*/
struct Inputs {
  Inputs();

  const Geometry *geometry;

  //! horizontal velocity components on the vertical grid (ghosted)
  const array::Array3D *u3;
  const array::Array3D *v3;
  //! vertical velocity relative to the bed
  const array::Array3D *w3;

  //! ice thickness change due to the surface mass balance during this step (m ice, positive
  //! = gain), as applied by GeometryEvolution
  const array::Scalar *top_surface_mass_balance;
  //! ice thickness change due to the basal mass balance during this step (m ice)
  const array::Scalar *bottom_surface_mass_balance;

  //! may be null
  const array::Scalar1 *no_model_mask;

  //! Throw if `geometry` is missing.
  void check() const;
  //! Throw if any of the fields needed by a transport model is missing.
  void check_transport() const;
};

//! A purely virtual class defining the interface of a PISM Debris Model.
class DebrisModel : public Component {
public:
  DebrisModel(std::shared_ptr<const Grid> g);
  DebrisModel(std::shared_ptr<const Grid> g, std::shared_ptr<DebrisModel> input);
  virtual ~DebrisModel() = default;

  void init(const Geometry &geometry);

  void update(const Inputs &inputs, double t, double dt);

  //! @brief Sets result to the mean debris, in "m".
  const array::Scalar& debris() const;


  void begin_pointwise_access() const;
  void end_pointwise_access() const;
  void init_timeseries(const std::vector<double> &ts) const;
  //! \brief Sets a pre-allocated N-element array "result" to the time-series
  //! of near-surface air temperature (kelvin) at the point i,j on the
  //! grid. Times (in years) are specified in ts. NB! Has to be surrounded by
  //! begin_pointwise_access() and end_pointwise_access()
  void debris_time_series(int i, int j, std::vector<double> &result) const;

protected:
  virtual void init_impl(const Geometry &geometry) = 0;
  virtual void update_impl(const Inputs &inputs, double t, double dt) = 0;

  virtual std::set<VariableMetadata> state_impl() const;

  virtual void write_state_impl(const OutputFile &output) const;

  virtual MaxTimestep max_timestep_impl(double t, const CFLData *cfl_data) const;

  virtual const array::Scalar& debris_impl() const;

  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;
  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual void debris_time_series_impl(int i, int j, std::vector<double> &result) const;

  virtual DiagnosticList spatial_diagnostics_impl() const;
  virtual TSDiagnosticList scalar_diagnostics_impl() const;

  mutable std::vector<double> m_ts_times;

  std::shared_ptr<DebrisModel> m_input_model;

  static std::shared_ptr<array::Scalar> allocate_debris(std::shared_ptr<const Grid> grid);
};

} // end of namespace debris
} // end of namespace pism

#endif  // PISM_DEBRIS_MODEL_HH
