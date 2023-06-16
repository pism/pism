// Copyright (C) 2008-2021 PISM Authors
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

#ifndef __PISMOceanModel_hh
#define __PISMOceanModel_hh

#include <memory>

#include "pism/util/Component.hh"

namespace pism {

//! @brief Ocean models and modifiers: provide sea level elevation,
//! melange back pressure, shelf base mass flux and shelf base
//! temperature.
namespace ocean {
//! A very rudimentary PISM ocean model.
class OceanModel : public Component {
public:
  // "modifier" constructor
  OceanModel(std::shared_ptr<const Grid> g, std::shared_ptr<OceanModel> input);
  // "model" constructor
  OceanModel(std::shared_ptr<const Grid> g);

  virtual ~OceanModel() = default;

  void init(const Geometry &geometry);

  void update(const Geometry &geometry, double t, double dt);

  const array::Scalar& shelf_base_temperature() const;
  const array::Scalar& shelf_base_mass_flux() const;
  const array::Scalar& average_water_column_pressure() const;

protected:
  virtual void init_impl(const Geometry &geometry);
  // provides default (pass-through) implementations for "modifiers"
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  virtual const array::Scalar& shelf_base_temperature_impl() const;
  virtual const array::Scalar& shelf_base_mass_flux_impl() const;
  virtual const array::Scalar& average_water_column_pressure_impl() const;

protected:
  std::shared_ptr<OceanModel> m_input_model;
  array::Scalar::Ptr m_water_column_pressure;

  static array::Scalar::Ptr allocate_shelf_base_temperature(std::shared_ptr<const Grid> g);
  static array::Scalar::Ptr allocate_shelf_base_mass_flux(std::shared_ptr<const Grid> g);
  static array::Scalar::Ptr allocate_water_column_pressure(std::shared_ptr<const Grid> g);

  static void compute_average_water_column_pressure(const Geometry &geometry,
                                                       double ice_density,
                                                       double water_density,
                                                       double g,
                                                       array::Scalar &result);

};

} // end of namespace ocean
} // end of namespace pism

#endif  // __PISMOceanModel_hh
