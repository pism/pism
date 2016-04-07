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

#include <cmath>                // for erf() in method 1 in putTempAtDepth()
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

namespace pism {

//! Read file and use heuristics to initialize PISM from typical 2d data available through remote sensing.
/*!
This procedure is called by the base class when option `-bootstrap` is used.

See chapter 4 of the User's Manual.  We read only 2D information from the bootstrap file.
 */
void IceModel::bootstrapFromFile(const std::string &filename) {

  // Bootstrap 2D fields:
  bootstrap_2d(filename);

  // Regrid 2D fields:
  regrid(2);

  // Check the consistency of geometry fields:
  updateSurfaceElevationAndMask();

  // save age, temperature, and enthalpy "revision numbers". If they
  // changed, then the corresponding field was initialized using
  // regridding.
  int age_revision = m_ice_age.get_state_counter(),
    temperature_revision = m_ice_temperature.get_state_counter(),
    enthalpy_revision = m_ice_enthalpy.get_state_counter();

  regrid(3);

  m_log->message(2,
             "bootstrapping 3D variables...\n");

  // Fill 3D fields using heuristics:
  {
    // set the initial age of the ice if appropriate
    if (m_config->get_boolean("do_age")) {
      if (age_revision == m_ice_age.get_state_counter()) {
        m_log->message(2,
                   " - setting initial age to %.4f years\n",
                   m_config->get_double("initial_age_of_ice_years"));
        m_ice_age.set(m_config->get_double("initial_age_of_ice_years", "seconds"));

      } else {
        m_log->message(2,
                   " - age of the ice was set already\n");
      }
    }

    if (m_config->get_boolean("do_cold_ice_methods") == true) {
      if (temperature_revision == m_ice_temperature.get_state_counter()) {
        m_log->message(2,
                   "getting surface B.C. from couplers...\n");
        init_step_couplers();

        // this call will set ice temperature
        putTempAtDepth();
      } else {
        m_log->message(2,
                   " - ice temperature was set already\n");
      }

      // make sure that enthalpy gets initialized:
      compute_enthalpy_cold(m_ice_temperature, m_ice_enthalpy);
      m_log->message(2,
                 " - ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");

    } else {
      // enthalpy mode
      if (enthalpy_revision == m_ice_enthalpy.get_state_counter()) {
        m_log->message(2,
                   "getting surface B.C. from couplers...\n");
        init_step_couplers();

        // this call will set ice enthalpy
        putTempAtDepth();
      } else {
        m_log->message(2,
                   " - ice enthalpy was set already\n");
      }
    }
  } // end of heuristics

  m_log->message(2, "done reading %s; bootstrapping done\n",filename.c_str());
}

void IceModel::bootstrap_2d(const std::string &filename) {

  PIO nc(m_grid->com, "guess_mode");
  nc.open(filename, PISM_READONLY);

  m_log->message(2,
             "bootstrapping by PISM default method from file %s\n", filename.c_str());

  // report on resulting computational box, rescale grid
  m_log->message(2,
             "  rescaling computational box for ice from -i file and\n"
             "    user options to dimensions:\n"
             "    [%6.2f km, %6.2f km] x [%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
             (m_grid->x0() - m_grid->Lx())/1000.0,
             (m_grid->x0() + m_grid->Lx())/1000.0,
             (m_grid->y0() - m_grid->Ly())/1000.0,
             (m_grid->y0() + m_grid->Ly())/1000.0,
             m_grid->Lz());

  std::string usurf_name;
  bool usurf_found = false, mask_found = false, usurf_found_by_std_name = false;
  nc.inq_var("usurf", "surface_altitude",
             usurf_found, usurf_name, usurf_found_by_std_name);
  mask_found = nc.inq_var("mask");

  std::string lon_name, lat_name;
  bool lon_found = false, lat_found = false,
    lon_found_by_std_name = false, lat_found_by_std_name = false;
  nc.inq_var("lon", "longitude", lon_found, lon_name, lon_found_by_std_name);
  nc.inq_var("lat", "latitude",  lat_found, lat_name, lat_found_by_std_name);

  nc.close();

  // now work through all the 2d variables, regridding if present and otherwise
  // setting to default values appropriately

  if (mask_found) {
    m_log->message(2,
               "  WARNING: 'mask' found; IGNORING IT!\n");
  }

  if (usurf_found) {
    m_log->message(2,
               "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n");
  }

  m_log->message(2,
             "  reading 2D model state variables by regridding ...\n");

  vLongitude.regrid(filename, OPTIONAL);
  if (not lon_found) {
    vLongitude.metadata().set_string("missing_at_bootstrap","true");
  }

  vLatitude.regrid(filename, OPTIONAL);
  if (not lat_found) {
    vLatitude.metadata().set_string("missing_at_bootstrap","true");
  }

  m_basal_melt_rate.regrid(filename, OPTIONAL,
                         m_config->get_double("bootstrapping_bmelt_value_no_var"));
  m_geothermal_flux.regrid(filename, OPTIONAL,
                         m_config->get_double("bootstrapping_geothermal_flux_value_no_var"));

  m_ice_thickness.regrid(filename, OPTIONAL,
                       m_config->get_double("bootstrapping_H_value_no_var"));
  // check the range of the ice thickness
  {
    Range thk_range = m_ice_thickness.range();

    if (thk_range.max >= m_grid->Lz() + 1e-6) {
      throw RuntimeError::formatted("Maximum ice thickness (%f meters)\n"
                                    "exceeds the height of the computational domain (%f meters).",
                                    thk_range.max, m_grid->Lz());
    }
  }

  if (m_config->get_boolean("part_grid")) {
    // Read the Href field from an input file. This field is
    // grid-dependent, so interpolating it from one grid to a
    // different one does not make sense in general.
    // (IceModel::Href_cleanup() will take care of the side effects of
    // such interpolation, though.)
    //
    // On the other hand, we need to read it in to be able to re-start
    // from a PISM output file using the -bootstrap option.
    vHref.regrid(filename, OPTIONAL, 0.0);
  }

  if (m_config->get_string("calving_methods").find("eigen_calving") != std::string::npos) {
    m_strain_rates.set(0.0);
  }

  if (m_config->get_boolean("ssa_dirichlet_bc")) {
    // Do not use Dirichlet B.C. anywhere if bc_mask is not present.
    m_ssa_dirichlet_bc_mask.regrid(filename, OPTIONAL, 0.0);
    // In the absence of u_ssa_bc and v_ssa_bc in the file the only B.C. that
    // makes sense is the zero Dirichlet B.C.
    m_ssa_dirichlet_bc_values.regrid(filename, OPTIONAL,  0.0);
  }

  // check if Lz is valid
  Range thk_range = m_ice_thickness.range();

  if (thk_range.max > m_grid->Lz()) {
    throw RuntimeError::formatted("Max. ice thickness (%3.3f m)\n"
                                  "exceeds the height of the computational domain (%3.3f m).",
                                  thk_range.max, m_grid->Lz());
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
`config.get_double("bootstrapping_temperature_heuristic"))`. Allowed values are `"smb"` and
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
void IceModel::putTempAtDepth() {

  m_log->message(2, " - filling ice temperatures using surface temps (and %s)\n",
             (m_config->get_string("bootstrapping_temperature_heuristic") == "quartic_guess"
              ? "quartic guess sans smb" : "mass balance for velocity estimate"));

  const bool do_cold = m_config->get_boolean("do_cold_ice_methods"),
             usesmb  = m_config->get_string("bootstrapping_temperature_heuristic") == "smb";
  const double
    ice_k = m_config->get_double("ice_thermal_conductivity"),
    melting_point_temp = m_config->get_double("water_melting_point_temperature"),
    ice_density = m_config->get_double("ice_density"),
    beta_CC_grad = m_config->get_double("beta_CC") * ice_density * m_config->get_double("standard_gravity"),
    KK = ice_k / (ice_density * m_config->get_double("ice_specific_heat_capacity"));

  assert(m_surface != NULL);

  // Note that surface->update was called in bootstrapFromFile()
  {
    m_surface->ice_surface_temperature(m_ice_surface_temp);
    if (usesmb == true) {
      m_surface->ice_surface_mass_flux(m_climatic_mass_balance);
      // convert from [kg m-2 s-1] to [m second-1]
      m_climatic_mass_balance.scale(1.0 / m_config->get_double("ice_density"));
    }
  }

  IceModelVec3 *result;
  if (do_cold) {
    result = &m_ice_temperature;
  } else {
    result = &m_ice_enthalpy;
  }

  IceModelVec::AccessList list;
  list.add(m_ice_surface_temp);
  list.add(m_climatic_mass_balance);
  list.add(m_ice_thickness);
  list.add(m_geothermal_flux);
  list.add(*result);

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double HH = m_ice_thickness(i,j),
        Ts = m_ice_surface_temp(i,j),
        gg = m_geothermal_flux(i,j);
      const unsigned int ks = m_grid->kBelowHeight(HH);

      double *T = result->get_column(i, j);

      // within ice
      if (usesmb == true) { // method 1:  includes surface mass balance in estimate
        const double mm = m_climatic_mass_balance(i,j);
        if (mm <= 0.0) { // negative or zero surface mass balance case: linear
          for (unsigned int k = 0; k < ks; k++) {
            const double z = m_grid->z(k),
              Tpmp = melting_point_temp - beta_CC_grad * (HH - z);
            T[k] = gg / ice_k * (HH - z) + Ts;
            T[k] = std::min(Tpmp,T[k]);
          }
        } else { // positive surface mass balance case
          const double C0 = (gg * sqrt(M_PI * HH * KK)) / (ice_k * sqrt(2.0 * mm)),
            gamma0 = sqrt(mm * HH / (2.0 * KK));

          for (unsigned int k = 0; k < ks; k++) {
            const double z = m_grid->z(k),
              Tpmp = melting_point_temp - beta_CC_grad * (HH - z);
            T[k] = Ts + C0 * (erf(gamma0) - erf(gamma0 * z / HH));
            T[k] = std::min(Tpmp,T[k]);
          }
        }
      } else { // method 2:  does not use smb
        const double beta = (4.0/21.0) * (gg / (2.0 * ice_k * HH * HH * HH)),
          alpha = (gg / (2.0 * HH * ice_k)) - 2.0 * HH * HH * beta;
        for (unsigned int k = 0; k < ks; k++) {
          const double depth = HH - m_grid->z(k),
            Tpmp = melting_point_temp - beta_CC_grad * depth,
            d2 = depth * depth;
          T[k] = std::min(Tpmp, Ts + alpha * d2 + beta * d2 * d2);
        }
      }

      // above ice
      for (unsigned int k = ks; k < m_grid->Mz(); k++) {
        T[k] = Ts;
      }

      // convert to enthalpy if that's what we are calculating
      if (not do_cold) {
        for (unsigned int k = 0; k < m_grid->Mz(); ++k) {
          const double depth = HH - m_grid->z(k);
          const double pressure = EC->pressure(depth);
          // reuse T to store enthalpy; assume that the ice is cold
          T[k]= EC->enthalpy_permissive(T[k], 0.0, pressure);
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->update_ghosts();
}


} // end of namespace pism
