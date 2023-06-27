// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2023 PISM Authors
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

#ifndef PISM_ICEBIN_SURFACE_MODEL_H
#define PISM_ICEBIN_SURFACE_MODEL_H

#include "pism/coupler/SurfaceModel.hh"

#include <memory>

namespace pism {
namespace icebin {

class IBSurfaceModel : public pism::surface::SurfaceModel {
public:
  IBSurfaceModel(std::shared_ptr<const pism::Grid> grid);

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double my_t, double my_dt);

  const array::Scalar& accumulation_impl() const;
  const array::Scalar& liquid_water_fraction_impl() const;
  const array::Scalar& mass_flux_impl() const;
  const array::Scalar& melt_impl() const;
  const array::Scalar& runoff_impl() const;
  const array::Scalar& temperature_impl() const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  MaxTimestep max_timestep_impl(double t) const;

  std::string m_input_file;

public:
  // ------ See icebin/contracts/modele_pism.cpp
  // Inputs from IceBin

  // Mass of ice being transferred GCM --> Ice Model
  pism::array::Scalar massxfer; // [kg m-2 s-1]
  // Enthalpy of ice being transferred Stieglitz --> Icebin
  pism::array::Scalar enthxfer; // [J m-2 s-1]

  // GCM's idea of energy transfer into ice sheet.
  // Used to compute mass/energy budget
  pism::array::Scalar deltah;

  // Temperature of the Dirichlet B.C.
  pism::array::Scalar ice_top_bc_temp;
  // Water content of the Dirichlet B.C.
  pism::array::Scalar ice_top_bc_wc;
};

} // end of namespace icebin
} // end of namespace pism

#endif /* PISM_ICEBIN_SURFACE_MODEL_H */
