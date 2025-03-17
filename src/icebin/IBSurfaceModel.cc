// Copyright (C) 2008-2016, 2023, 2024, 2025 PISM Authors
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

#include "pism/util/MaxTimestep.hh"
#include "pism/util/io/File.hh"
#include "pism/icebin/IBSurfaceModel.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace icebin {

IBSurfaceModel::IBSurfaceModel(std::shared_ptr<const pism::Grid> grid)
    : SurfaceModel(grid),
      massxfer(grid, "massxfer"),
      enthxfer(grid, "enthxfer"),
      deltah(grid, "deltah"),
      ice_top_bc_temp(grid, "ice_top_bc_temp"),
      ice_top_bc_wc(grid, "ice_top_bc_wc") {

  massxfer.metadata(0)
      .long_name("Mass of ice being transferred Stieglitz --> Icebin")
      .units("kg m^-2 s^-1")
      .standard_name("land_ice_surface_specific_mass_balance");

  enthxfer.metadata(0)
      .long_name("Enthalpy of ice being transferred Stieglitz --> Icebin")
      .units("W m^-2");

  // ------- Used only for mass/energy budget
  deltah.metadata(0)
      .long_name(
          "enthalpy of constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate")
      .units("W m^-2");

  // ------- Dirichlet Bondary condition derived from deltah
  ice_top_bc_temp.metadata(0).long_name("Temperature of the Dirichlet B.C.").units("kelvin");
  ice_top_bc_wc.metadata(0).long_name("Water content of the Dirichlet B.C.").units("1");
}

void IBSurfaceModel::init_impl(const Geometry &geometry) {
  (void)geometry;

  m_log->message(2, "* Initializing the IceBin interface surface model IBSurfaceModel.\n"
                    "  IceBin changes its state when surface conditions change.\n");

  for (auto *v : {&massxfer, &enthxfer, &deltah, &ice_top_bc_temp, &ice_top_bc_wc}) {
    v->set(0.0);
  }
}

MaxTimestep IBSurfaceModel::max_timestep_impl(double t) const {
  (void)t;
  return {};
}

void IBSurfaceModel::update_impl(const Geometry &geometry, double t, double dt) {
  (void)geometry;
  (void)t;
  (void)dt;

  // compute naive estimates of accumulation, melt, and runoff
  dummy_accumulation(massxfer, *m_accumulation);
  dummy_melt(massxfer, *m_melt);
  dummy_runoff(massxfer, *m_runoff);
}

const array::Scalar &IBSurfaceModel::liquid_water_fraction_impl() const {
  return ice_top_bc_wc;
}

const array::Scalar &IBSurfaceModel::mass_flux_impl() const {
  return massxfer;
}

const array::Scalar &IBSurfaceModel::temperature_impl() const {
  return ice_top_bc_temp;
}

const array::Scalar& IBSurfaceModel::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar& IBSurfaceModel::melt_impl() const {
  return *m_melt;
}

const array::Scalar& IBSurfaceModel::runoff_impl() const {
  return *m_runoff;
}

void IBSurfaceModel::define_model_state_impl(const File &output) const {
  for (const auto *v : {&massxfer, &enthxfer, &deltah, &ice_top_bc_temp, &ice_top_bc_wc}) {
    v->define(output, io::PISM_DOUBLE);
  }
}

void IBSurfaceModel::write_model_state_impl(const File &output) const {
  for (const auto *v : {&massxfer, &enthxfer, &deltah, &ice_top_bc_temp, &ice_top_bc_wc}) {
    v->write(output);
  }
}

} // namespace icebin
} // namespace pism
