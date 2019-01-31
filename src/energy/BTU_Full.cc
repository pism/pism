/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
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
#include "BTU_Full.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"
#include "BedrockColumn.hh"

namespace pism {
namespace energy {


BTU_Full::BTU_Full(IceGrid::ConstPtr g, const BTUGrid &grid)
  : BedThermalUnit(g),
    m_bootstrapping_needed(false) {

  m_k = m_config->get_double("energy.bedrock_thermal_conductivity");

  const double
    rho = m_config->get_double("energy.bedrock_thermal_density"),
    c   = m_config->get_double("energy.bedrock_thermal_specific_heat_capacity");
  // build constant diffusivity for heat equation
  m_D   = m_k / (rho * c);

  // validate Lbz
  if (grid.Lbz <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Invalid bedrock thermal layer depth: %f m",
                                  grid.Lbz);
  }

  // validate Mbz
  if (grid.Mbz < 2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Invalid number of layers of the bedrock thermal layer: %d",
                                  grid.Mbz);
  }

  {
    m_Mbz = grid.Mbz;
    m_Lbz = grid.Lbz;

    std::map<std::string, std::string> attrs;
    attrs["units"]     = "m";
    attrs["long_name"] = "Z-coordinate in bedrock";
    attrs["axis"]      = "Z";
    attrs["positive"]  = "up";

    std::vector<double> z(m_Mbz);
    double dz = m_Lbz / (m_Mbz - 1);
    for (unsigned int k = 0; k < m_Mbz; ++k) {
      z[k] = -m_Lbz + k * dz;
    }
    z.back() = 0.0;
    m_temp.create(m_grid, "litho_temp", "zb", z, attrs);

    m_temp.set_attrs("model_state",
                    "lithosphere (bedrock) temperature, in BTU_Full",
                    "K", "");
    m_temp.metadata().set_double("valid_min", 0.0);
  }

  m_column.reset(new BedrockColumn("bedrock_column", *m_config, vertical_spacing(), Mz()));
}

BTU_Full::~BTU_Full() {
  // empty
}


//! \brief Initialize the bedrock thermal unit.
void BTU_Full::init_impl(const InputOptions &opts) {

  m_log->message(2, "* Initializing the bedrock thermal unit...\n");

  // 2D initialization. Takes care of the flux through the bottom surface of the thermal layer.
  BedThermalUnit::init_impl(opts);

  // Initialize the temperature field.
  {
    // store the current "revision number" of the temperature field
    const int temp_revision = m_temp.state_counter();

    if (opts.type == INIT_RESTART) {
      PIO input_file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);

      if (input_file.inq_var("litho_temp")) {
        m_temp.read(input_file, opts.record);
      }
      // otherwise the bedrock temperature is either interpolated from a -regrid_file or filled
      // using bootstrapping (below)
    }

    regrid("bedrock thermal layer", m_temp, REGRID_WITHOUT_REGRID_VARS);

    if (m_temp.state_counter() == temp_revision) {
      m_bootstrapping_needed = true;
    } else {
      m_bootstrapping_needed = false;
    }
  }

  update_flux_through_top_surface();
}


/** Returns the vertical spacing used by the bedrock grid.
 */
double BTU_Full::vertical_spacing_impl() const {
  return m_Lbz / (m_Mbz - 1.0);
}

unsigned int BTU_Full::Mz_impl() const {
  return m_Mbz;
}


double BTU_Full::depth_impl() const {
  return m_Lbz;
}

void BTU_Full::define_model_state_impl(const PIO &output) const {
  m_bottom_surface_flux.define(output);
  m_temp.define(output);
}

void BTU_Full::write_model_state_impl(const PIO &output) const {
  m_bottom_surface_flux.write(output);
  m_temp.write(output);
}

MaxTimestep BTU_Full::max_timestep_impl(double t) const {
  (void) t;
  // No time step restriction: we are using an unconditionally stable method.
  return MaxTimestep("bedrock thermal layer");
}


/** Perform a step of the bedrock thermal model.
*/
void BTU_Full::update_impl(const IceModelVec2S &bedrock_top_temperature,
                           double t, double dt) {
  (void) t;

  if (m_bootstrapping_needed) {
    bootstrap(bedrock_top_temperature);
    m_bootstrapping_needed = false;
  }

  if (dt < 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dt < 0 is not allowed");
  }

  IceModelVec::AccessList list{&m_temp, &m_bottom_surface_flux, &bedrock_top_temperature};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double *T = m_temp.get_column(i, j);

      m_column->solve(dt, m_bottom_surface_flux(i, j), bedrock_top_temperature(i, j),
                      T,  // input
                      T); // output

      // Check that T is positive:
      for (unsigned int k = 0; k < m_Mbz; ++k) {
        if (T[k] <= 0.0) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "invalid bedrock temperature: %f Kelvin at %d,%d,%d",
                                        T[k], i, j, k);
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  update_flux_through_top_surface();
}

/*! Computes the heat flux from the bedrock thermal layer upward into the
  ice/bedrock interface:
  \f[G_0 = -k_b \frac{\partial T_b}{\partial z}\big|_{z=0}.\f]
  Uses the second-order finite difference expression
  \f[\frac{\partial T_b}{\partial z}\big|_{z=0} \approx \frac{3 T_b(0) - 4 T_b(-\Delta z) + T_b(-2\Delta z)}{2 \Delta z}\f]
  where \f$\Delta z\f$ is the equal spacing in the bedrock.

  The above expression only makes sense when `Mbz` = `temp.n_levels` >= 3.
  When `Mbz` = 2 we use first-order differencing.  When temp was not created,
  the `Mbz` <= 1 cases, we return the stored geothermal flux.
*/
void BTU_Full::update_flux_through_top_surface() {

  if (m_bootstrapping_needed) {
    m_top_surface_flux.copy_from(m_bottom_surface_flux);
    return;
  }

  double dz = this->vertical_spacing();
  const int k0  = m_Mbz - 1;  // Tb[k0] = ice/bed interface temp, at z=0

  IceModelVec::AccessList list{&m_temp, &m_top_surface_flux};

  if (m_Mbz >= 3) {

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Tb = m_temp.get_column(i, j);
      m_top_surface_flux(i, j) = - m_k * (3 * Tb[k0] - 4 * Tb[k0-1] + Tb[k0-2]) / (2 * dz);
    }

  } else {

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Tb = m_temp.get_column(i, j);
      m_top_surface_flux(i, j) = - m_k * (Tb[k0] - Tb[k0-1]) / dz;
    }

  }
}

const IceModelVec3Custom& BTU_Full::temperature() const {
  if (m_bootstrapping_needed) {
    throw RuntimeError(PISM_ERROR_LOCATION, "bedrock temperature is not available (bootstrapping is needed)");
  }

  return m_temp;
}

void BTU_Full::bootstrap(const IceModelVec2S &bedrock_top_temperature) {

  m_log->message(2,
                "  bootstrapping to fill lithosphere temperatures in the bedrock thermal layer\n"
                "  using temperature at the top bedrock surface and geothermal flux\n"
                "  (bedrock temperature is linear in depth)...\n");

  double dz = this->vertical_spacing();
  const int k0 = m_Mbz - 1; // Tb[k0] = ice/bedrock interface temp

  IceModelVec::AccessList list{&bedrock_top_temperature, &m_bottom_surface_flux, &m_temp};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *Tb = m_temp.get_column(i, j); // Tb points into temp memory

    Tb[k0] = bedrock_top_temperature(i, j);
    for (int k = k0-1; k >= 0; k--) {
      Tb[k] = Tb[k+1] + dz * m_bottom_surface_flux(i, j) / m_k;
    }
  }

  m_temp.inc_state_counter();     // mark as modified
}

} // end of namespace energy
} // end of namespace pism
