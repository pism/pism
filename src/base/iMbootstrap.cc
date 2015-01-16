// Copyright (C) 2004-2015 Jed Brown, Nathan Shemonski, Ed Bueler and
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

#include "iceModel.hh"
#include "PIO.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include <cmath>                // for erf() in method 1 in putTempAtDepth()
#include <cassert>
#include "error_handling.hh"

namespace pism {

//! Read file and use heuristics to initialize PISM from typical 2d data available through remote sensing.
/*! 
This procedure is called by the base class when option `-boot_file` is used.

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
  int age_revision = age3.get_state_counter(),
    temperature_revision = T3.get_state_counter(),
    enthalpy_revision = Enth3.get_state_counter();

  regrid(3);

  verbPrintf(2, grid.com,
             "bootstrapping 3D variables...\n");

  // Fill 3D fields using heuristics:
  {
    // set the initial age of the ice if appropriate
    if (config.get_flag("do_age")) {
      if (age_revision == age3.get_state_counter()) {
        verbPrintf(2, grid.com,
                   " - setting initial age to %.4f years\n",
                   config.get("initial_age_of_ice_years"));
        age3.set(config.get("initial_age_of_ice_years", "years", "seconds"));

      } else {
        verbPrintf(2, grid.com,
                   " - age of the ice was set already\n");
      }
    }

    if (config.get_flag("do_cold_ice_methods") == true) {
      if (temperature_revision == T3.get_state_counter()) {
        verbPrintf(2, grid.com,
                   "getting surface B.C. from couplers...\n");
        init_step_couplers();

        // this call will set ice temperature
        putTempAtDepth();
      } else {
        verbPrintf(2, grid.com,
                   " - ice temperature was set already\n");
      }

      // make sure that enthalpy gets initialized:
      compute_enthalpy_cold(T3, Enth3);
      verbPrintf(2, grid.com,
                 " - ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");

    } else {
      // enthalpy mode
      if (enthalpy_revision == Enth3.get_state_counter()) {
        verbPrintf(2, grid.com,
                   "getting surface B.C. from couplers...\n");
        init_step_couplers();

        // this call will set ice enthalpy
        putTempAtDepth();
      } else {
        verbPrintf(2, grid.com,
                   " - ice enthalpy was set already\n");
      }
    }
  } // end of heuristics

  verbPrintf(2, grid.com, "done reading %s; bootstrapping done\n",filename.c_str());
}

void IceModel::bootstrap_2d(const std::string &filename) {

  PIO nc(grid, "guess_mode");
  nc.open(filename, PISM_READONLY);

  verbPrintf(2, grid.com, 
             "bootstrapping by PISM default method from file %s\n", filename.c_str());

  // report on resulting computational box, rescale grid
  verbPrintf(2, grid.com, 
             "  rescaling computational box for ice from -boot_file file and\n"
             "    user options to dimensions:\n"
             "    [%6.2f km, %6.2f km] x [%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
             (grid.x0() - grid.Lx())/1000.0,
             (grid.x0() + grid.Lx())/1000.0,
             (grid.y0() - grid.Ly())/1000.0,
             (grid.y0() + grid.Ly())/1000.0,
             grid.Lz());

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
    verbPrintf(2, grid.com, 
               "  WARNING: 'mask' found; IGNORING IT!\n");
  }

  if (usurf_found) {
    verbPrintf(2, grid.com, 
               "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n");
  }

  verbPrintf(2, grid.com, 
             "  reading 2D model state variables by regridding ...\n");

  vLongitude.regrid(filename, OPTIONAL);
  if (not lon_found) {
    vLongitude.metadata().set_string("missing_at_bootstrap","true");
  }

  vLatitude.regrid(filename, OPTIONAL);
  if (not lat_found) {
    vLatitude.metadata().set_string("missing_at_bootstrap","true");
  }

  basal_melt_rate.regrid(filename, OPTIONAL,
                         config.get("bootstrapping_bmelt_value_no_var"));
  geothermal_flux.regrid(filename, OPTIONAL,
                         config.get("bootstrapping_geothermal_flux_value_no_var"));

  ice_thickness.regrid(filename, OPTIONAL,
                       config.get("bootstrapping_H_value_no_var"));
  // check the range of the ice thickness
  {
    Range thk_range = ice_thickness.range();

    if (thk_range.max >= grid.Lz() + 1e-6) {
      throw RuntimeError::formatted("Maximum ice thickness (%f meters)\n"
                                    "exceeds the height of the computational domain (%f meters).",
                                    thk_range.max, grid.Lz());
    }
  }

  if (config.get_flag("part_grid")) {
    // Read the Href field from an input file. This field is
    // grid-dependent, so interpolating it from one grid to a
    // different one does not make sense in general.
    // (IceModel::Href_cleanup() will take care of the side effects of
    // such interpolation, though.)
    //
    // On the other hand, we need to read it in to be able to re-start
    // from a PISM output file using the -boot_file option.
    vHref.regrid(filename, OPTIONAL, 0.0);
  }

  if (config.get_string("calving_methods").find("eigen_calving") != std::string::npos) {
    strain_rates.set(0.0);
  }

  if (config.get_flag("ssa_dirichlet_bc")) {
    // Do not use Dirichlet B.C. anywhere if bcflag is not present.
    vBCMask.regrid(filename, OPTIONAL, 0.0);
    // In the absence of u_ssa_bc and v_ssa_bc in the file the only B.C. that
    // makes sense is the zero Dirichlet B.C.
    vBCvel.regrid(filename, OPTIONAL,  0.0);
  }

  // check if Lz is valid
  Range thk_range = ice_thickness.range();

  if (thk_range.max > grid.Lz()) {
    throw RuntimeError::formatted("Max. ice thickness (%3.3f m)\n"
                                  "exceeds the height of the computational domain (%3.3f m).",
                                  thk_range.max, grid.Lz());
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
`config.get("bootstrapping_temperature_heuristic"))`. Allowed values are `"smb"` and `"quartic_guess"`.

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

  verbPrintf(2, grid.com, " - filling ice temperatures using surface temps (and %s)\n",
             (config.get_string("bootstrapping_temperature_heuristic") == "quartic_guess"
              ? "quartic guess sans smb" : "mass balance for velocity estimate"));

  const bool do_cold = config.get_flag("do_cold_ice_methods"),
             usesmb  = config.get_string("bootstrapping_temperature_heuristic") == "smb";
  const double
    ice_k = config.get("ice_thermal_conductivity"),
    melting_point_temp = config.get("water_melting_point_temperature"),
    ice_density = config.get("ice_density"),
    beta_CC_grad = config.get("beta_CC") * ice_density * config.get("standard_gravity"),
    KK = ice_k / (ice_density * config.get("ice_specific_heat_capacity"));

  assert(surface != NULL);

  // Note that surface->update was called in bootstrapFromFile()
  {
    surface->ice_surface_temperature(ice_surface_temp);
    if (usesmb == true) {
      surface->ice_surface_mass_flux(climatic_mass_balance);
      // convert from [kg m-2 s-1] to [m / s]
      climatic_mass_balance.scale(1.0 / config.get("ice_density"));
    }
  }

  IceModelVec3 *result;
  if (do_cold) {
    result = &T3;
  } else {
    result = &Enth3;
  }

  IceModelVec::AccessList list;
  list.add(ice_surface_temp);
  list.add(climatic_mass_balance);
  list.add(ice_thickness);
  list.add(geothermal_flux);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double HH = ice_thickness(i,j),
      Ts = ice_surface_temp(i,j),
      gg = geothermal_flux(i,j);
    const unsigned int ks = grid.kBelowHeight(HH);

    double *T = NULL;
    result->getInternalColumn(i, j, &T);

    // within ice
    if (usesmb == true) { // method 1:  includes surface mass balance in estimate
      const double mm = climatic_mass_balance(i,j);
      if (mm <= 0.0) { // negative or zero surface mass balance case: linear
        for (unsigned int k = 0; k < ks; k++) {
          const double z = grid.z(k),
            Tpmp = melting_point_temp - beta_CC_grad * (HH - z);
          T[k] = gg / ice_k * (HH - z) + Ts;
          T[k] = std::min(Tpmp,T[k]);
        }
      } else { // positive surface mass balance case
        const double C0 = (gg * sqrt(M_PI * HH * KK)) / (ice_k * sqrt(2.0 * mm)),
          gamma0 = sqrt(mm * HH / (2.0 * KK));

        for (unsigned int k = 0; k < ks; k++) {
          const double z = grid.z(k),
            Tpmp = melting_point_temp - beta_CC_grad * (HH - z);
          T[k] = Ts + C0 * (erf(gamma0) - erf(gamma0 * z / HH));
          T[k] = std::min(Tpmp,T[k]);
        }
      }
    } else { // method 2:  does not use smb
      const double beta = (4.0/21.0) * (gg / (2.0 * ice_k * HH * HH * HH)),
        alpha = (gg / (2.0 * HH * ice_k)) - 2.0 * HH * HH * beta;
      for (unsigned int k = 0; k < ks; k++) {
        const double depth = HH - grid.z(k),
          Tpmp = melting_point_temp - beta_CC_grad * depth,
          d2 = depth * depth;
        T[k] = std::min(Tpmp, Ts + alpha * d2 + beta * d2 * d2);
      }
    }

    // above ice
    for (unsigned int k = ks; k < grid.Mz(); k++) {
      T[k] = Ts;
    }

    // convert to enthalpy if that's what we are calculating
    if (do_cold == false) {
      for (unsigned int k = 0; k < grid.Mz(); ++k) {
        const double depth = HH - grid.z(k);
        const double pressure = EC->getPressureFromDepth(depth);
        // reuse T to store enthalpy; assume that the ice is cold
        T[k]= EC->getEnthPermissive(T[k], 0.0, pressure);
      }
    }

  }


  result->update_ghosts();
}


} // end of namespace pism
