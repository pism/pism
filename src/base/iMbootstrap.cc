// Copyright (C) 2004-2016 Jed Brown, Nathan Shemonski, Ed Bueler and
// Constantine Khroulev
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

#include <cmath>                // for erf() in method 1 in bootstrap_ice_temperature()
#include <cassert>
#include <gsl/gsl_math.h>       // M_PI

#include "iceModel.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "enthalpyConverter.hh"
#include "base/energy/BedThermalUnit.hh"
#include "base/energy/utilities.hh"

namespace pism {

void IceModel::bootstrap_2d(const PIO &input_file) {

  m_log->message(2, "bootstrapping from file '%s'...\n", input_file.inq_filename().c_str());

  std::string usurf_name;
  bool usurf_found = false, mask_found = false, usurf_found_by_std_name = false;
  input_file.inq_var("usurf", "surface_altitude",
                     usurf_found, usurf_name, usurf_found_by_std_name);
  mask_found = input_file.inq_var("mask");

  std::string lon_name, lat_name;
  bool lon_found = false, lat_found = false,
    lon_found_by_std_name = false, lat_found_by_std_name = false;
  input_file.inq_var("lon", "longitude", lon_found, lon_name, lon_found_by_std_name);
  input_file.inq_var("lat", "latitude",  lat_found, lat_name, lat_found_by_std_name);

  // now work through all the 2d variables, regridding if present and otherwise
  // setting to default values appropriately

  if (mask_found) {
    m_log->message(2, "  WARNING: 'mask' found; IGNORING IT!\n");
  }

  if (usurf_found) {
    m_log->message(2, "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n");
  }

  m_log->message(2, "  reading 2D model state variables by regridding ...\n");

  m_longitude.regrid(input_file, OPTIONAL);
  if (not lon_found) {
    m_longitude.metadata().set_string("missing_at_bootstrap","true");
  }

  m_latitude.regrid(input_file, OPTIONAL);
  if (not lat_found) {
    m_latitude.metadata().set_string("missing_at_bootstrap","true");
  }

  m_ice_thickness.regrid(input_file, OPTIONAL,
                         m_config->get_double("bootstrapping.defaults.ice_thickness"));
  // check the range of the ice thickness
  {
    Range thk_range = m_ice_thickness.range();

    if (thk_range.max >= m_grid->Lz() + 1e-6) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Maximum ice thickness (%f meters)\n"
                                    "exceeds the height of the computational domain (%f meters).",
                                    thk_range.max, m_grid->Lz());
    }
  }

  if (m_config->get_boolean("geometry.part_grid.enabled")) {
    // Read the Href field from an input file. This field is
    // grid-dependent, so interpolating it from one grid to a
    // different one does not make sense in general.
    // (IceModel::Href_cleanup() will take care of the side effects of
    // such interpolation, though.)
    //
    // On the other hand, we need to read it in to be able to re-start
    // from a PISM output file using the -bootstrap option.
    m_Href.regrid(input_file, OPTIONAL, 0.0);
  }

  if (m_config->get_boolean("stress_balance.ssa.dirichlet_bc")) {
    // Do not use Dirichlet B.C. anywhere if bc_mask is not present.
    m_ssa_dirichlet_bc_mask.regrid(input_file, OPTIONAL, 0.0);
    // In the absence of u_ssa_bc and v_ssa_bc in the file the only B.C. that
    // makes sense is the zero Dirichlet B.C.
    m_ssa_dirichlet_bc_values.regrid(input_file, OPTIONAL,  0.0);
  }

  // check if Lz is valid
  Range thk_range = m_ice_thickness.range();

  if (thk_range.max > m_grid->Lz()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Max. ice thickness (%3.3f m)\n"
                                  "exceeds the height of the computational domain (%3.3f m).",
                                  thk_range.max, m_grid->Lz());
  }
}

} // end of namespace pism
