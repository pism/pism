// Copyright (C) 2023 PISM Authors
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

#ifndef PISM_PYOCEANMODEL_H
#define PISM_PYOCEANMODEL_H

#include "pism/coupler/ocean/CompleteOceanModel.hh"
#include <memory>

namespace pism {
namespace ocean {

/*!
 * The base class for ocean models that are implemented in Python
 */
class PyOceanModel {
public:
  virtual ~PyOceanModel() = default;

  std::shared_ptr<pism::array::Scalar> shelf_base_temperature;
  std::shared_ptr<array::Scalar> shelf_base_mass_flux;
  std::shared_ptr<array::Scalar> water_column_pressure;

  /*!
   * Allocate data members. We need this to be able to test a Python implementation of an
   * ocean model *without* the rest of PISM.
   */
  void allocate(std::shared_ptr<const Grid> grid);

  /*!
   * Maximum time step the model can take at time `t` (in seconds)
   */
  virtual MaxTimestep max_timestep(double t) const;

  /*!
   * Initialize the model state from formulas, by reading from an input file, etc.
   */
  virtual void init(const Geometry &geometry);

  /*!
   * Update the state of the model by taking a time step from `t` to `t + dt` (in seconds).
   *
   * Assumes that the time step length `dt` is allowed at the time `t`.
   */
  virtual void update(const Geometry &geometry, double t, double dt);

  /*!
   * Define model state variables and set their attributes
   */
  virtual void define_model_state(const File &output) const;

  /*!
   * Write model state variables and set their attributes
   */
  virtual void write_model_state(const File &output) const;
};

//! The adapter class for Python ocean models
/*!
 * We need this class because SWIG cannot create a wrapper for the OceanModel class that
 * can be used as a base class for Python classes. (Specifically: it does not support
 * methods that return a reference, e.g. `const array::Scalar& shelf_base_mass_flux_impl()`.)
 */
class PyOceanModelAdapter : public CompleteOceanModel {
public:
  PyOceanModelAdapter(std::shared_ptr<const Grid> grid, std::shared_ptr<PyOceanModel> implementation);
  virtual ~PyOceanModelAdapter() = default;

private:
  MaxTimestep max_timestep_impl(double t) const;
  void update_impl(const Geometry &geometry, double my_t, double my_dt);
  void init_impl(const Geometry &geometry);

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  std::shared_ptr<PyOceanModel> m_impl;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* PISM_PYOCEANMODEL_H */
