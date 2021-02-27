// Copyright (C) 2008-2020 PISM Authors
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

#include <gsl/gsl_math.h>       // GSL_NAN

#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Vars.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/icebin/IBSurfaceModel.hh"

namespace pism {
namespace icebin {

///// Constant-in-time surface model for accumulation,
///// ice surface temperature parameterized as in PISM-IBSurfaceModel dependent on latitude and surface elevation


IBSurfaceModel::IBSurfaceModel(IceGrid::ConstPtr g)
  : SurfaceModel(g),
    icebin_wflux(m_grid, "icebin_wflux", WITHOUT_GHOSTS),
    icebin_deltah(m_grid, "icebin_deltah", WITHOUT_GHOSTS),
    icebin_massxfer(m_grid, "icebin_massxfer", WITHOUT_GHOSTS),
    icebin_enthxfer(m_grid, "icebin_enthxfer", WITHOUT_GHOSTS),
    surface_temp(m_grid, "surface_temp", WITHOUT_GHOSTS)
{

  printf("BEGIN IBSurfaceModel::allocate_IBSurfaceModel()\n");
  icebin_wflux.set_attrs("climate_state",
                         "constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
                         "kg m-2 s-1", "kg m-2 year-1", "land_ice_surface_specific_mass_balance", 0);

  icebin_deltah.set_attrs(
      "climate_state",
      "enthalpy of constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
      "W m-2", "W m-2", "", 0);

  icebin_massxfer.set_attrs(
      "climate_state",
      "enthalpy of constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
      "kg m-2 s-1", "kg m-2 s-1", "", 0);

  icebin_enthxfer.set_attrs("climate_state", "constant-in-time heat flux through top surface",
                            "W m-2", "W m-2", "", 0);

  // This variable is computed from the inputs above.
  surface_temp.set_attrs("climate_state", "Temperature to use for Dirichlet B.C. at surface",
                         "K", "K", "", 0);

  printf("END IBSurfaceModel::allocate_IBSurfaceModel()\n");
}

void IBSurfaceModel::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the IceBin interface surface model IBSurfaceModel.\n"
                    "  IceBin changes its state when surface conditions change.\n");

  // find PISM input file to read data from:
  m_input_file = process_input_options(m_grid->com, m_config).filename;

  // It doesn't matter what we set this to, it will be re-set later.
  icebin_wflux.set(0.0);
  icebin_deltah.set(0.0);
  icebin_massxfer.set(0.0);
  icebin_enthxfer.set(0.0);
  surface_temp.set(0.0);

  _initialized = true;
}

MaxTimestep IBSurfaceModel::max_timestep_impl(double t) const {
  (void)t;
  return MaxTimestep("surface icebin");
}

void IBSurfaceModel::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) t;
  (void) dt;
}

const IceModelVec2S &IBSurfaceModel::mass_flux_impl() const {
  return icebin_massxfer;
}

const IceModelVec2S &IBSurfaceModel::temperature_impl() const {
  return surface_temp;
}

void IBSurfaceModel::define_model_state_impl(const File &output) const {
  SurfaceModel::define_model_state_impl(output);
  icebin_enthxfer.define(output);
  icebin_wflux.define(output);
  icebin_deltah.define(output);
  icebin_massxfer.define(output);
  surface_temp.define(output);
}

void IBSurfaceModel::write_model_state_impl(const File &output) const {
  SurfaceModel::write_model_state_impl(output);
  icebin_enthxfer.write(output);
  icebin_wflux.write(output);
  icebin_deltah.write(output);
  icebin_massxfer.write(output);
  surface_temp.write(output);
}

} // end of namespace surface
} // end of namespace pism
