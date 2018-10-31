/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#include "utilities.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/Logger.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/pism_utilities.hh"
#include "bootstrapping.hh"

namespace pism {
namespace energy {

//! Compute ice enthalpy from temperature temperature by assuming the ice has zero liquid fraction.
/*!
First this method makes sure the temperatures is at most the pressure-melting
value, before computing the enthalpy for that temperature, using zero liquid
fraction.

Because of how EnthalpyConverter::pressure() works, the energy
content in the air is set to the value that ice would have if it a chunk of it
occupied the air; the atmosphere actually has much lower energy content.  It is
done this way for regularity (i.e. dEnth/dz computations).
*/
void compute_enthalpy_cold(const IceModelVec3 &temperature,
                           const IceModelVec2S &ice_thickness,
                           IceModelVec3 &result) {

  IceGrid::ConstPtr grid = result.grid();
  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  IceModelVec::AccessList list{&temperature, &result, &ice_thickness};

  const unsigned int Mz = grid->Mz();
  const std::vector<double> &z = grid->z();

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *Tij = temperature.get_column(i,j);
    double *Enthij = result.get_column(i,j);

    for (unsigned int k = 0; k < Mz; ++k) {
      const double depth = ice_thickness(i, j) - z[k]; // FIXME issue #15
      Enthij[k] = EC->enthalpy_permissive(Tij[k], 0.0, EC->pressure(depth));
    }
  }

  result.inc_state_counter();

  result.update_ghosts();
}

void compute_temperature(const IceModelVec3 &enthalpy,
                         const IceModelVec2S &ice_thickness,
                         IceModelVec3 &result) {

  IceGrid::ConstPtr grid = result.grid();
  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  IceModelVec::AccessList list{&enthalpy, &ice_thickness, &result};

  const unsigned int Mz = grid->Mz();
  const std::vector<double> &z = grid->z();

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      *E = enthalpy.get_column(i, j),
      H  = ice_thickness(i, j);
    double *T = result.get_column(i, j);

    for (unsigned int k = 0; k < Mz; ++k) {
      const double depth = H - z[k]; // FIXME issue #15
      T[k] = EC->temperature(E[k], EC->pressure(depth));
    }
  }

  result.inc_state_counter();

  result.update_ghosts();
}

//! Compute `result` (enthalpy) from `temperature` and liquid fraction.
void compute_enthalpy(const IceModelVec3 &temperature,
                      const IceModelVec3 &liquid_water_fraction,
                      const IceModelVec2S &ice_thickness,
                      IceModelVec3 &result) {

  IceGrid::ConstPtr grid = result.grid();
  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  IceModelVec::AccessList list{&temperature, &liquid_water_fraction, &ice_thickness, &result};

  const unsigned int Mz = grid->Mz();
  const std::vector<double> &z = grid->z();

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *T     = temperature.get_column(i,j);
    const double *omega = liquid_water_fraction.get_column(i,j);
    double       *E     = result.get_column(i,j);

    for (unsigned int k = 0; k < Mz; ++k) {
      const double depth = ice_thickness(i,j) - z[k]; // FIXME issue #15
      E[k] = EC->enthalpy_permissive(T[k], omega[k], EC->pressure(depth));
    }
  }

  result.update_ghosts();

  result.inc_state_counter();
}

//! Compute the liquid fraction corresponding to enthalpy and ice_thickness.
void compute_liquid_water_fraction(const IceModelVec3 &enthalpy,
                                   const IceModelVec2S &ice_thickness,
                                   IceModelVec3 &result) {

  IceGrid::ConstPtr grid = result.grid();

  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  result.set_name("liqfrac");
  result.metadata(0).set_name("liqfrac");
  result.set_attrs("diagnostic",
                   "liquid water fraction in ice (between 0 and 1)",
                   "1", "", 0);

  IceModelVec::AccessList list{&result, &enthalpy, &ice_thickness};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Enthij = enthalpy.get_column(i,j);
      double *omegaij = result.get_column(i,j);

      for (unsigned int k=0; k < grid->Mz(); ++k) {
        const double depth = ice_thickness(i,j) - grid->z(k); // FIXME issue #15
        omegaij[k] = EC->water_fraction(Enthij[k],EC->pressure(depth));
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.inc_state_counter();
}

//! @brief Compute the CTS field, CTS = E/E_s(p), from `ice_enthalpy` and `ice_thickness`, and put
//! in `result`.
/*!
 * The actual cold-temperate transition surface (CTS) is the level set CTS = 1.
 *
 * Does not communicate ghosts for IceModelVec3 result.
 */
void compute_cts(const IceModelVec3 &ice_enthalpy,
                 const IceModelVec2S &ice_thickness,
                 IceModelVec3 &result) {

  IceGrid::ConstPtr grid = result.grid();
  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  result.set_name("cts");
  result.metadata(0).set_name("cts");
  result.set_attrs("diagnostic",
                   "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
                   "", "", 0);

  IceModelVec::AccessList list{&ice_enthalpy, &ice_thickness, &result};

  const unsigned int Mz = grid->Mz();
  const std::vector<double> &z = grid->z();

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *CTS  = result.get_column(i,j);
    const double *enthalpy = ice_enthalpy.get_column(i,j);

    for (unsigned int k = 0; k < Mz; ++k) {
      const double depth = ice_thickness(i,j) - z[k]; // FIXME issue #15
      CTS[k] = enthalpy[k] / EC->enthalpy_cts(EC->pressure(depth));
    }
  }

  result.inc_state_counter();
}

//! Computes the total ice enthalpy in J.
/*!
  Units of the specific enthalpy field are J kg-1.  We integrate
  \f$E(t,x,y,z)\f$ over the entire ice fluid region \f$\Omega(t)\f$, multiplying
  by the density to get units of energy:
  \f[ E_{\text{total}}(t) = \int_{\Omega(t)} E(t,x,y,z) \rho_i \,dx\,dy\,dz. \f]
*/
double total_ice_enthalpy(double thickness_threshold,
                          const IceModelVec3 &ice_enthalpy,
                          const IceModelVec2S &ice_thickness) {
  double enthalpy_sum = 0.0;

  IceGrid::ConstPtr grid = ice_enthalpy.grid();
  Config::ConstPtr config = grid->ctx()->config();

  auto cell_area = grid->cell_area();

  const std::vector<double> &z = grid->z();

  IceModelVec::AccessList list{&ice_enthalpy, &ice_thickness};
  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = ice_thickness(i, j);

      if (H >= thickness_threshold) {
        const int ks = grid->kBelowHeight(H);

        const double
          *E   = ice_enthalpy.get_column(i, j);

        for (int k = 0; k < ks; ++k) {
          enthalpy_sum += cell_area * E[k] * (z[k+1] - z[k]);
        }
        enthalpy_sum += cell_area * E[ks] * (H - z[ks]);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  enthalpy_sum *= config->get_double("constants.ice.density");

  return GlobalSum(grid->com, enthalpy_sum);
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

  IceGrid::ConstPtr      grid   = result.grid();
  Context::ConstPtr      ctx    = grid->ctx();
  Config::ConstPtr       config = ctx->config();
  Logger::ConstPtr       log    = ctx->log();
  EnthalpyConverter::Ptr EC     = ctx->enthalpy_converter();

  const bool use_smb  = config->get_string("bootstrapping.temperature_heuristic") == "smb";

  if (use_smb) {
    log->message(2,
                 " - Filling 3D ice temperatures using surface temperature"
                 " (and mass balance for velocity estimate)...\n");

  } else {
    log->message(2,
                 " - Filling 3D ice temperatures using surface temperature"
                 " (and a quartic guess without SMB)...\n");
  }

  const double
    ice_k       = config->get_double("constants.ice.thermal_conductivity"),
    ice_density = config->get_double("constants.ice.density"),
    ice_c       = config->get_double("constants.ice.specific_heat_capacity"),
    K           = ice_k / (ice_density * ice_c),
    T_min       = config->get_double("energy.minimum_allowed_temperature"),
    T_melting   = config->get_double("constants.fresh_water.melting_point_temperature",
                                     "Kelvin");

  IceModelVec::AccessList list{&ice_surface_temp, &surface_mass_balance,
      &ice_thickness, &basal_heat_flux, &result};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        T_surface = std::min(ice_surface_temp(i, j), T_melting),
        H         = ice_thickness(i, j),
        G         = basal_heat_flux(i, j);

      if (G < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "geothermal flux G(%d,%d) = %f < 0.0 %s",
                                      i, j, G,
                                      basal_heat_flux.metadata().get_string("units").c_str());
      }

      if (T_surface < T_min) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "T_surface(%d,%d) = %f < T_min = %f Kelvin",
                                      i, j, T_surface, T_min);
      }

      const unsigned int ks = grid->kBelowHeight(H);

      double *T = result.get_column(i, j);

      // within ice
      if (use_smb) { // method 1:  includes surface mass balance in estimate

        // Convert SMB from "kg m-2 s-1" to "m second-1".
        const double SMB = surface_mass_balance(i, j) / ice_density;

        for (unsigned int k = 0; k < ks; k++) {
          const double z = grid->z(k);
          T[k] = ice_temperature_guess_smb(EC, H, z, T_surface, G, ice_k, K, SMB);
        }

      } else { // method 2: a quartic guess; does not use SMB

        for (unsigned int k = 0; k < ks; k++) {
          const double z = grid->z(k);
          T[k] = ice_temperature_guess(EC, H, z, T_surface, G, ice_k);
        }

      }

      // Make sure that resulting temperatures are not too low.
      for (unsigned int k = 0; k < ks; k++) {
        if (T[k] < T_min) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "T(%d,%d,%d) = %f < T_min = %f Kelvin",
                                        i, j, k, T[k], T_min);
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

} // end of namespace energy
} // end of namespace pism
