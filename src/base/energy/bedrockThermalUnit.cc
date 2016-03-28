// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Ed Bueler and Constantine Khroulev
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

#include <gsl/gsl_math.h>

#include "bedrockThermalUnit.hh"
#include "base/util/io/PIO.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include <cassert>
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace energy {

BedThermalUnit::BedThermalUnit(IceGrid::ConstPtr g)
    : Component_TS(g) {

  const Logger &log = *g->ctx()->log();

  m_upward_flux.create(m_grid, "bheatflx", WITHOUT_GHOSTS);
  // PROPOSED standard_name = lithosphere_upward_heat_flux
  m_upward_flux.set_attrs("diagnostic", "upward geothermal flux at bedrock surface",
                          "W m-2", "");
  m_upward_flux.metadata().set_string("glaciological_units", "mW m-2");
  m_upward_flux.write_in_glaciological_units = true;

  // build constant diffusivity for heat equation
  m_bed_rho = m_config->get_double("bedrock_thermal_density");
  m_bed_c   = m_config->get_double("bedrock_thermal_specific_heat_capacity");
  m_bed_k   = m_config->get_double("bedrock_thermal_conductivity");
  m_bed_D   = m_bed_k / (m_bed_rho * m_bed_c);

  m_Mbz = (int)m_config->get_double("grid_Mbz");
  m_Lbz = (int)m_config->get_double("grid_Lbz");
  m_input_file.clear();

  // FIXME: Move the code processing command-line options elsewhere,
  // possibly making Mbz and Lbz arguments of the constructor. It's
  // good to validate Lbz and Mbz here, though.
  {
    options::String input_file("-i", "PISM input file name");
    bool boot = options::Bool("-bootstrap", "enable bootstrapping heuristics");

    options::Integer Mbz_option("-Mbz", "number of levels in bedrock thermal layer", m_Mbz);
    m_Mbz = Mbz_option;

    options::Real Lbz("-Lbz",
                      "depth (thickness) of bedrock thermal layer, in meters", m_Lbz);
    m_Lbz = Lbz;

    if (input_file.is_set() and not boot) {
      m_input_file = input_file;
      options::ignored(log, "-Mbz");
      options::ignored(log, "-Lbz");

      // If we're initializing from a file we need to get the number of bedrock
      // levels and the depth of the bed thermal layer from it:
      PIO nc(m_grid->com, "guess_mode");

      nc.open(m_input_file, PISM_READONLY);

      bool exists = nc.inq_var("litho_temp");

      if (exists) {
        grid_info info(nc, "litho_temp", m_sys, m_grid->periodicity());

        m_Mbz = info.z_len;
        m_Lbz = -info.z_min;
      } else {
        // override values we got using config.get_double() in the constructor
        m_Mbz = 1;
        m_Lbz = 0;
      }

      nc.close();
    } else {
      // Bootstrapping

      if (Mbz_option.is_set() && m_Mbz == 1) {
        options::ignored(log, "-Lbz");
        m_Lbz = 0;
      } else if (Mbz_option.is_set() ^ Lbz.is_set()) {
        throw RuntimeError("please specify both -Mbz and -Lbz");
      }
    }

    // actual allocation

    // validate Lbz and Mbz:
    if ((m_Lbz <= 0.0) && (m_Mbz > 1)) {
      throw RuntimeError("BedThermalUnit can not be created with negative or zero Lbz value\n"
                         "and more than one layers");
    }

    if (m_Mbz > 1) {
      std::map<std::string, std::string> attrs;
      attrs["units"] = "m";
      attrs["long_name"] = "Z-coordinate in bedrock";
      attrs["axis"] = "Z";
      attrs["positive"] = "up";

      std::vector<double> z(m_Mbz);
      double dz = m_Lbz / (m_Mbz - 1);
      for (unsigned int k = 0; k < m_Mbz; ++k) {
        z[k] = -m_Lbz + k * dz;
      }
      z.back() = 0;
      m_temp.create(m_grid, "litho_temp", "zb", z, attrs);

      m_temp.set_attrs("model_state",
                     "lithosphere (bedrock) temperature, in BedThermalUnit",
                     "K", "");
      m_temp.metadata().set_double("valid_min", 0.0);
    }
  }
}

BedThermalUnit::~BedThermalUnit() {
  // empty
}


//! \brief Initialize the bedrock thermal unit.
void BedThermalUnit::init(bool &bootstrapping_needed) {
  grid_info g;

  // first assume that we don't need to bootstrap
  bootstrapping_needed = false;

  // store the current "revision number" of the temperature field
  int temp_revision = m_temp.get_state_counter();

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the bedrock thermal unit... setting constants...\n");

  // If we're using a minimal model, then we're done:
  if (!m_temp.was_created()) {
    m_log->message(2,
               "  minimal model for lithosphere: stored geothermal flux applied to ice base ...\n");
    return;
  }

  if (not m_input_file.empty()) {
    PIO nc(m_grid->com, "guess_mode");

    nc.open(m_input_file, PISM_READONLY);
    bool exists = nc.inq_var("litho_temp");

    if (exists) {
      const unsigned int last_record = nc.inq_nrecords("litho_temp", "", m_sys) - 1;
      m_temp.read(m_input_file, last_record);
    }

    nc.close();
  }

  if (m_temp.was_created() == true) {
    regrid("BedThermalUnit", m_temp, REGRID_WITHOUT_REGRID_VARS);
  }

  if (m_temp.get_state_counter() == temp_revision) {
    bootstrapping_needed = true;
  }
}

/** Returns the vertical spacing used by the bedrock grid.
 *
 * Special case: returns 0 if the bedrock thermal layer has thickness
 * zero.
 */
double BedThermalUnit::vertical_spacing() const {
  if (m_temp.was_created() == true) {
    return m_Lbz / (m_Mbz - 1.0);
  } else {
    return 0.0;
  }
}

unsigned int BedThermalUnit::Mbz() const {
  return m_Mbz;
}

void BedThermalUnit::add_vars_to_output_impl(const std::string &/*keyword*/, std::set<std::string> &result) {
  if (m_temp.was_created()) {
    result.insert(m_temp.metadata().get_string("short_name"));
  }
}

void BedThermalUnit::define_variables_impl(const std::set<std::string> &vars,
                                                const PIO &nc, IO_Type nctype) {
  if (m_temp.was_created()) {
    if (set_contains(vars, m_temp.metadata().get_string("short_name"))) {
      m_temp.define(nc, nctype);
    }
  }
}

void BedThermalUnit::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  if (m_temp.was_created()) {
    if (set_contains(vars, m_temp.metadata().get_string("short_name"))) {
      m_temp.write(nc);
    }
  }
}

void BedThermalUnit::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                          std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
  dict["hfgeoubed"] = Diagnostic::Ptr(new BTU_geothermal_flux_at_ground_level(this));
  (void)ts_dict;
}


/*! Because the grid for the bedrock thermal layer is equally-spaced, and because
the heat equation being solved in the bedrock is time-invariant (%e.g. no advection
at evolving velocity and no time-dependence to physical constants), the explicit
time-stepping can compute the maximum stable time step easily.  The basic scheme
is
        \f[T_k^{n+1} = T_k^n + R (T_{k-1}^n - 2 T_k^n + T_{k+1}^n)\f]
where
        \f[R = \frac{k \Delta t}{\rho c \Delta z^2} = \frac{D \Delta t}{\Delta z^2}.\f]
The stability condition is that the coefficients of temperatures on the right are
all nonnegative, equivalently \f$1-2R\ge 0\f$ or \f$R\le 1/2\f$ or
        \f[\Delta t \le \frac{\Delta z^2}{2 D}.\f]
This is a formula for the maximum stable timestep.  For more, see [\ref MortonMayers].

The above describes the general case where Mbz > 1.
 */
MaxTimestep BedThermalUnit::max_timestep_impl(double t) {
  (void) t;

  if (m_temp.was_created()) {
    double dzb = this->vertical_spacing();
    // max dt from stability; in seconds
    return MaxTimestep(dzb * dzb / (2.0 * m_bed_D));
  } else {
    return MaxTimestep();
  }
}


/** Perform a step of the bedrock thermal model.

@todo The old scheme had better stability properties, as follows:

Because there is no advection, the simplest centered implicit (backward Euler) scheme is easily "bombproof" without choosing \f$\lambda\f$, or other complications.  It has this scaled form,
\anchor bedrockeqn
\f[ -R_b T_{k-1}^{n+1} + \left(1 + 2 R_b\right) T_k^{n+1} - R_b T_{k+1}^{n+1}
         = T_k^n, \tag{bedrockeqn} \f]
where
  \f[ R_b = \frac{k_b \Delta t}{\rho_b c_b \Delta z^2}. \f]
This is unconditionally stable for a pure bedrock problem, and has a maximum principle, without any further qualification [\ref MortonMayers].

@todo Now a trapezoid rule could be used
*/
void BedThermalUnit::update_impl(double my_t, double my_dt) {

  if (not m_temp.was_created()) {
    update_upward_geothermal_flux();
    return;  // in this case we are up to date
  }

  // as a derived class of Component_TS, has t,dt members which keep track
  // of last update time-interval; so we do some checks ...
  // CHECK: has the desired time-interval already been dealt with?
  if ((fabs(my_t - m_t) < 1e-12) && (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  // CHECK: is the desired time interval a forward step?; backward heat equation not good!
  if (my_dt < 0) {
     throw RuntimeError("BedThermalUnit::update() does not allow negative timesteps");
  }

  // CHECK: is desired time-interval equal to [my_t,my_t+my_dt] where my_t = t + dt?
  if ((!gsl_isnan(m_t)) && (!gsl_isnan(m_dt))) { // this check should not fire on first use
    bool contiguous = true;

    if (fabs(m_t + m_dt) < 1) {
      if (fabs(my_t - (m_t + m_dt)) >= 1e-12) { // check if the absolute difference is small
        contiguous = false;
      }
    } else {
      if (fabs(my_t - (m_t + m_dt)) / (m_t + m_dt) >= 1e-12) { // check if the relative difference is small
        contiguous = false;
      }
    }

    if (not contiguous) {
      throw RuntimeError::formatted("BedThermalUnit::update() requires next update to be contiguous with last;\n"
                                    "  stored:     t = %f s,    dt = %f s\n"
                                    "  desired: my_t = %f s, my_dt = %f s",
                                    m_t,m_dt,my_t,my_dt); }
  }

  // CHECK: is desired time-step too long?
  MaxTimestep my_max_dt = max_timestep(my_t);
  if (my_max_dt.is_finite() and my_max_dt.value() < my_dt) {
     throw RuntimeError("BedThermalUnit::update() thinks you asked for too big a timestep.");
  }

  // o.k., we have checked; we are going to do the desired timestep!
  m_t  = my_t;
  m_dt = my_dt;

  // Get pointers to fields owned by IceModel.
  const IceModelVec2S
    *bedtoptemp = m_grid->variables().get_2d_scalar("bedtoptemp"),
    *ghf        = m_grid->variables().get_2d_scalar("bheatflx");

  double dzb = this->vertical_spacing();
  const int  k0  = m_Mbz - 1;          // Tb[k0] = ice/bed interface temp, at z=0

  const double bed_R  = m_bed_D * my_dt / (dzb * dzb);

  std::vector<double> Tbnew(m_Mbz);

  IceModelVec::AccessList list;
  list.add(m_temp);
  list.add(*ghf);
  list.add(*bedtoptemp);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *Tbold = m_temp.get_column(i,j); // Tbold actually points into temp memory
    Tbold[k0] = (*bedtoptemp)(i,j);  // sets Dirichlet explicit-in-time b.c. at top of bedrock column

    const double Tbold_negone = Tbold[1] + 2 * (*ghf)(i,j) * dzb / m_bed_k;
    Tbnew[0] = Tbold[0] + bed_R * (Tbold_negone - 2 * Tbold[0] + Tbold[1]);
    for (int k = 1; k < k0; k++) { // working upward from base
      Tbnew[k] = Tbold[k] + bed_R * (Tbold[k-1] - 2 * Tbold[k] + Tbold[k+1]);
    }
    Tbnew[k0] = (*bedtoptemp)(i,j);

    m_temp.set_column(i,j,&Tbnew[0]); // copy from Tbnew into temp memory
  }

  update_upward_geothermal_flux();
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
void BedThermalUnit::update_upward_geothermal_flux() {

  if (not m_temp.was_created()) {
    m_upward_flux.copy_from(*m_grid->variables().get_2d_scalar("bheatflx"));
    return;
  }

  double dzb = this->vertical_spacing();
  const int k0  = m_Mbz - 1;  // Tb[k0] = ice/bed interface temp, at z=0

  IceModelVec::AccessList list;
  list.add(m_temp);
  list.add(m_upward_flux);

  if (m_Mbz >= 3) {

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Tb = m_temp.get_column(i,j);
      m_upward_flux(i,j) = - m_bed_k * (3 * Tb[k0] - 4 * Tb[k0-1] + Tb[k0-2]) / (2 * dzb);
    }

  } else {

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Tb = m_temp.get_column(i,j);
      m_upward_flux(i,j) = - m_bed_k * (Tb[k0] - Tb[k0-1]) / dzb;
    }

  }
}

const IceModelVec2S& BedThermalUnit::upward_geothermal_flux() const {
  return m_upward_flux;
}

void BedThermalUnit::bootstrap() {

  if (m_Mbz < 2) {
    return;
  }

  m_log->message(2,
             "  bootstrapping to fill lithosphere temperatures in bedrock thermal layers,\n"
             "    using provided bedtoptemp and a linear function from provided geothermal flux ...\n");

  double dzb = this->vertical_spacing();
  const int k0 = m_Mbz-1; // Tb[k0] = ice/bedrock interface temp

  // Get pointers to fields owned by IceModel.
  const IceModelVec2S
    *bedtoptemp = m_grid->variables().get_2d_scalar("bedtoptemp"),
    *ghf        = m_grid->variables().get_2d_scalar("bheatflx");

  IceModelVec::AccessList list;
  list.add(*bedtoptemp);
  list.add(*ghf);
  list.add(m_temp);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *Tb = m_temp.get_column(i,j); // Tb points into temp memory
    Tb[k0] = (*bedtoptemp)(i,j);
    for (int k = k0-1; k >= 0; k--) {
      Tb[k] = Tb[k+1] + dzb * (*ghf)(i,j) / m_bed_k;
    }
  }

  m_temp.inc_state_counter();     // mark as modified
}

BTU_geothermal_flux_at_ground_level::BTU_geothermal_flux_at_ground_level(BedThermalUnit *m)
  : Diag<BedThermalUnit>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "hfgeoubed"));
  set_attrs("upward geothermal flux at ground (top of the bedrock) level",
            "",                 // no standard name
            "W m-2", "W m-2", 0);
}

IceModelVec::Ptr BTU_geothermal_flux_at_ground_level::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hfgeoubed", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  result->copy_from(model->upward_geothermal_flux());

  return result;
}

} // end of namespace energy
} // end of namespace pism
