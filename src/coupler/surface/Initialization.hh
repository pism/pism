/* Copyright (C) 2016, 2017, 2018, 2019, 2022 PISM Authors
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

#ifndef PSINITIALIZATION_H
#define PSINITIALIZATION_H

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

/*! Surface model "modifier" that helps with initialization.
 *
 * This modifier saves *all* fields a surface model provides as a part of the model state and
 * re-loads them during initialization so that they are available *before* the first time step in a
 * re-started run.
 *
 * It is
 *
 * - not visible to the user,
 * - is added automatically, and
 * - does not have a corresponding "keyword" in surface::Factory.
 */
class InitializationHelper : public SurfaceModel {
public:
  InitializationHelper(std::shared_ptr<const IceGrid> g, std::shared_ptr<SurfaceModel> in);
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const array::Scalar &layer_mass_impl() const;
  const array::Scalar &liquid_water_fraction_impl() const;
  const array::Scalar &temperature_impl() const;
  const array::Scalar &mass_flux_impl() const;
  const array::Scalar &layer_thickness_impl() const;

  const array::Scalar& accumulation_impl() const;
  const array::Scalar& melt_impl() const;
  const array::Scalar& runoff_impl() const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

private:
  // store pointers to fields so that we can iterate over them
  std::vector<array::Array*> m_variables;
  // storage
  array::Scalar m_mass_flux;
  array::Scalar m_temperature;
  // the rest of the field are inherited from SurfaceModel
};

} // end of namespace surface
} // end of namespace pism


#endif /* PSINITIALIZATION_H */
