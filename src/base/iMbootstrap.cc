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

  m_basal_melt_rate.regrid(input_file, OPTIONAL,
                           m_config->get_double("bootstrapping.defaults.bmelt"));

  m_ice_thickness.regrid(input_file, OPTIONAL,
                         m_config->get_double("bootstrapping.defaults.ice_thickness"));
  // check the range of the ice thickness
  {
    Range thk_range = m_ice_thickness.range();

    if (thk_range.max >= m_grid->Lz() + 1e-6) {
      throw RuntimeError::formatted("Maximum ice thickness (%f meters)\n"
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

  if (m_config->get_string("calving.methods").find("eigen_calving") != std::string::npos) {
    m_strain_rates.set(0.0);
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
    throw RuntimeError::formatted("Max. ice thickness (%3.3f m)\n"
                                  "exceeds the height of the computational domain (%3.3f m).",
                                  thk_range.max, m_grid->Lz());
  }
}

//! Fill 3D fields using heuristics.
void IceModel::bootstrap_3d() {

  // set the initial age of the ice if appropriate
  if (m_config->get_boolean("age.enabled")) {
    m_log->message(2, " - setting initial age to %.4f years\n",
                   m_config->get_double("age.initial_value"));

    m_ice_age.set(m_config->get_double("age.initial_value", "seconds"));
  }

  {
    m_surface->ice_surface_temperature(m_ice_surface_temp);
    m_surface->ice_surface_mass_flux(m_climatic_mass_balance);
  }

  if (m_config->get_boolean("energy.temperature_based")) {
    // set ice temperature:
    bootstrap_ice_temperature(m_ice_thickness,
                              m_ice_surface_temp,
                              m_climatic_mass_balance,
                              m_btu->flux_through_top_surface(),
                              m_ice_temperature);

    // use temperature to initialize enthalpy:
    compute_enthalpy_cold(m_ice_temperature, m_ice_thickness, m_ice_enthalpy);

    m_log->message(2, " - ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");
  } else {
    // enthalpy mode

    bootstrap_ice_enthalpy(m_ice_thickness,
                           m_ice_surface_temp,
                           m_climatic_mass_balance,
                           m_btu->flux_through_top_surface(),
                           m_ice_enthalpy);
  }
}

//! Create a temperature field within the ice from provided ice thickness, surface temperature, surface mass balance, and geothermal flux.
/*!
In bootstrapping we need to determine initial values for the temperature within
the ice (and the bedrock).  There are various data available at bootstrapping,
but not the 3D temperature field needed as initial values for the temperature.  Here
we take a "guess" based on an assumption of steady state and a simple model of
the vertical velocity in the column.  The rule is certainly heuristic but it
seems to work well anyway.

The result is *not* the temperature field which is in steady state with the ice
dynamics.  Spinup is most-definitely needed in many applications.  Such spinup
usually starts from the temperature field computed by this procedure and then
runs for a long time (e.g. \f$10^4\f$ to \f$10^6\f$ years), possibly with fixed
geometry, to get closer to thermomechanically-coupled equilibrium.

Consider a horizontal grid point.  Suppose the surface temperature
\f$T_s\f$, surface mass balance \f$m\f$, and geothermal flux \f$g\f$ are given at that location.
Within the column denote the temperature by \f$T(z)\f$ at height \f$z\f$ above
the base of the ice.  Suppose the column of ice has height \f$H\f$, the ice
thickness.

There are two alternative bootstrap methods determined by the configuration parameter
`config.get_double("bootstrapping.temperature_heuristic"))`. Allowed values are `"smb"` and
`"quartic_guess"`.

1. If the `smb` method is chosen, which is the default, and if \f$m>0\f$,
then the method sets the ice
temperature to the solution of the steady problem [\ref Paterson]
  \f[\rho_i c w \frac{\partial T}{\partial z} = k_i \frac{\partial^2 T}{\partial z^2} \qquad \text{with boundary conditions} \qquad T(H) = T_s \quad \text{and} \quad \frac{\partial T}{\partial z}(0) = - \frac{g}{k_i}, \f]
where the vertical velocity is linear between the surface value \f$w=-m\f$ and
a velocity of zero at the base:
  \f[w(z) = - m z / H.\f]
(Note that because \f$m>0\f$, this vertical velocity is downward.)
This is a two-point boundary value problem for a linear ODE.  In fact, if
\f$K = k_i / (\rho_i c)\f$ then we can write the ODE as
  \f[K T'' + \frac{m z}{H} T' = 0.\f]
Then let
  \f[C_0 = \frac{g \sqrt{\pi H K}}{k_i \sqrt{2 m}}, \qquad \gamma_0 = \sqrt{\frac{mH}{2K}}.\f]
(Note \f$\gamma_0\f$ is, up to a constant, the square root of the Peclet number
[\ref Paterson]; compare [\ref vanderWeletal2013].)  The solution to the
two-point boundary value problem is then
  \f[T(z) = T_s + C_0 \left(\operatorname{erf}(\gamma_0) - \operatorname{erf}\left(\gamma_0 \frac{z}{H}\right)\right).\f]
If `usesmb` is true and \f$m \le 0\f$, then the velocity in the column, relative
to the base, is taken to be zero.  Thus the solution is
  \f[ T(z) = \frac{g}{k_i} \left( H - z \right) + T_s, \f]
a straight line whose slope is determined by the geothermal flux and whose value
at the ice surface is the surface temperature, \f$T(H) = T_s\f$.
2. If the `quartic_guess` method is chosen, the "quartic guess" formula which was in older
versions of PISM is used.  Namely, within the ice we set
\f[T(z) = T_s + \alpha (H-z)^2 + \beta (H-z)^4\f]
where \f$\alpha,\beta\f$ are chosen so that
\f[\frac{\partial T}{\partial z}\Big|_{z=0} = - \frac{g}{k_i} \qquad \text{and} \qquad \frac{\partial T}{\partial z}\Big|_{z=H/4} = - \frac{g}{2 k_i}.\f]
The purpose of the second condition is that when ice is advecting downward then
the temperature gradient is much larger in roughly the bottom quarter of the
ice column.  However, without the surface mass balance, much less the solution
of the stress balance equations, we cannot estimate the vertical velocity, so
we make such a rough guess.

In either case the temperature within the ice is not allowed to exceed the
pressure-melting temperature.

We set \f$T(z)=T_s\f$ above the top of the ice.

This method determines \f$T(0)\f$, the ice temperature at the ice base.  This
temperature is used by BedThermalUnit::bootstrap() to determine a
bootstrap temperature profile in the bedrock.
*/
void bootstrap_ice_temperature(const IceModelVec2S &ice_thickness,
                               const IceModelVec2S &ice_surface_temp,
                               const IceModelVec2S &surface_mass_balance,
                               const IceModelVec2S &basal_heat_flux,
                               IceModelVec3 &result) {

  IceGrid::ConstPtr      grid   = result.get_grid();
  Context::ConstPtr      ctx    = grid->ctx();
  Config::ConstPtr       config = ctx->config();
  Logger::ConstPtr       log    = ctx->log();
  EnthalpyConverter::Ptr EC     = ctx->enthalpy_converter();

  const bool use_smb  = config->get_string("bootstrapping.temperature_heuristic") == "smb";

  if (use_smb) {
    log->message(2,
                 " - filling 3D ice temperatures using surface temperature"
                 " (and mass balance for velocity estimate)\n");

  } else {
    log->message(2,
                 " - filling 3D ice temperatures using surface temperature"
                 " (and a quartic guess without SMB)\n");
  }

  const double
    ice_k       = config->get_double("constants.ice.thermal_conductivity"),
    ice_density = config->get_double("constants.ice.density"),
    ice_c       = config->get_double("constants.ice.specific_heat_capacity"),
    K           = ice_k / (ice_density * ice_c);

  IceModelVec::AccessList list;
  list.add(ice_surface_temp);
  list.add(surface_mass_balance);
  list.add(ice_thickness);
  list.add(basal_heat_flux);
  list.add(result);

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        T_surface = ice_surface_temp(i, j),
        H         = ice_thickness(i, j),
        G         = basal_heat_flux(i, j);

      const unsigned int ks = grid->kBelowHeight(H);

      double *T = result.get_column(i, j);

      // within ice
      if (use_smb) { // method 1:  includes surface mass balance in estimate

        // Convert SMB from "kg m-2 s-1" to "m second-1".
        const double SMB = surface_mass_balance(i, j) / ice_density;

        if (SMB <= 0.0) {
          // negative or zero surface mass balance case: linear
          for (unsigned int k = 0; k < ks; k++) {
            const double
              depth = H - grid->z(k),
              Tpmp  = EC->melting_temperature(EC->pressure(depth));

            T[k] = G / ice_k * depth + T_surface;
            T[k] = std::min(Tpmp, T[k]);
          }
        } else {
          // positive surface mass balance case
          const double
            C0     = (G * sqrt(M_PI * H * K)) / (ice_k * sqrt(2.0 * SMB)),
            gamma0 = sqrt(SMB * H / (2.0 * K));

          for (unsigned int k = 0; k < ks; k++) {
            const double
              z    = grid->z(k),
              Tpmp = EC->melting_temperature(EC->pressure(H - z));

            T[k] = T_surface + C0 * (erf(gamma0) - erf(gamma0 * z / H));
            T[k] = std::min(Tpmp,T[k]);
          }
        }
      } else { // method 2:  does not use SMB
        const double
          beta = (4.0/21.0) * (G / (2.0 * ice_k * H * H * H)),
          alpha = (G / (2.0 * H * ice_k)) - 2.0 * H * H * beta;

        for (unsigned int k = 0; k < ks; k++) {
          const double
            depth = H - grid->z(k),
            Tpmp  = EC->melting_temperature(EC->pressure(depth)),
            d2    = depth * depth;

          T[k] = std::min(Tpmp, T_surface + alpha * d2 + beta * d2 * d2);
        }
      }

      // above ice
      for (unsigned int k = ks; k < grid->Mz(); k++) {
        T[k] = T_surface;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();
}

void bootstrap_ice_enthalpy(const IceModelVec2S &ice_thickness,
                            const IceModelVec2S &ice_surface_temp,
                            const IceModelVec2S &surface_mass_balance,
                            const IceModelVec2S &basal_heat_flux,
                            IceModelVec3 &result) {

  bootstrap_ice_temperature(ice_thickness, ice_surface_temp,
                            surface_mass_balance, basal_heat_flux,
                            result);

  compute_enthalpy_cold(result, ice_thickness, result);
}


} // end of namespace pism
