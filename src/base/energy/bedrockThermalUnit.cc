// Copyright (C) 2011, 2012, 2013, 2014 Ed Bueler and Constantine Khroulev
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

#include "bedrockThermalUnit.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include <assert.h>
#include "PISMConfig.hh"
#include "error_handling.hh"

namespace pism {

BedThermalUnit::BedThermalUnit(const IceGrid &g)
    : Component_TS(g) {
  bedtoptemp = NULL;
  ghf        = NULL;

  // build constant diffusivity for heat equation
  bed_rho = m_config.get("bedrock_thermal_density");
  bed_c   = m_config.get("bedrock_thermal_specific_heat_capacity");
  bed_k   = m_config.get("bedrock_thermal_conductivity");
  bed_D   = bed_k / (bed_rho * bed_c);

  m_Mbz = (int)m_config.get("grid_Mbz");
  m_Lbz = (int)m_config.get("grid_Lbz");
  m_input_file.clear();

  // FIXME: Move the code processing command-line options elsewhere,
  // possibly making Mbz and Lbz arguments of the constructor. It's
  // good to validate Lbz and Mbz here, though.
  {
    options::String i("-i", "PISM input file name");

    options::Integer Mbz("-Mbz", "number of levels in bedrock thermal layer", m_Mbz);
    m_Mbz = Mbz;

    options::Real Lbz("-Lbz",
                      "depth (thickness) of bedrock thermal layer, in meters", m_Lbz);
    m_Lbz = Lbz;

    if (i.is_set()) {
      m_input_file = i;
      ignore_option(m_grid.com, "-Mbz");
      ignore_option(m_grid.com, "-Lbz");

      // If we're initializing from a file we need to get the number of bedrock
      // levels and the depth of the bed thermal layer from it:
      PIO nc(m_grid, "guess_mode");

      nc.open(m_input_file, PISM_READONLY);

      bool exists = nc.inq_var("litho_temp");

      if (exists) {
        grid_info info = nc.inq_grid_info("litho_temp", m_grid.periodicity());

        m_Mbz = info.z_len;
        m_Lbz = -info.z_min;
      } else {
        // override values we got using config.get() in the constructor
        m_Mbz = 1;
        m_Lbz = 0;
      }

      nc.close();
    } else {
      // Bootstrapping

      if (Mbz.is_set() && m_Mbz == 1) {
        ignore_option(m_grid.com, "-Lbz");
        m_Lbz = 0;
      } else if (Mbz.is_set() ^ Lbz.is_set()) {
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
      temp.create(m_grid, "litho_temp", "zb", z, attrs);

      temp.set_attrs("model_state",
                     "lithosphere (bedrock) temperature, in BedThermalUnit",
                     "K", "");
      temp.metadata().set_double("valid_min", 0.0);
    }
  }
}


//! \brief Initialize the bedrock thermal unit.
void BedThermalUnit::init(bool &bootstrapping_needed) {
  grid_info g;

  // first assume that we don't need to bootstrap
  bootstrapping_needed = false;

  // store the current "revision number" of the temperature field
  int temp_revision = temp.get_state_counter();

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2,m_grid.com,
             "* Initializing the bedrock thermal unit... setting constants...\n");

  // Get pointers to fields owned by IceModel.
  bedtoptemp = m_grid.variables().get_2d_scalar("bedtoptemp");
  ghf = m_grid.variables().get_2d_scalar("bheatflx");

  // If we're using a minimal model, then we're done:
  if (!temp.was_created()) {
    verbPrintf(2,m_grid.com,
               "  minimal model for lithosphere: stored geothermal flux applied to ice base ...\n");
    return;
  }

  if (not m_input_file.empty()) {
    PIO nc(m_grid, "guess_mode");

    nc.open(m_input_file, PISM_READONLY);
    bool exists = nc.inq_var("litho_temp");

    if (exists) {
      const unsigned int last_record = nc.inq_nrecords("litho_temp", "") - 1;
      temp.read(m_input_file, last_record);
    }

    nc.close();
  }

  if (temp.was_created() == true) {
    regrid("BedThermalUnit", &temp, REGRID_WITHOUT_REGRID_VARS);
  }

  if (temp.get_state_counter() == temp_revision) {
    bootstrapping_needed = true;
  }
}

/** Returns the vertical spacing used by the bedrock grid.
 *
 * Special case: returns 0 if the bedrock thermal layer has thickness
 * zero.
 */
double BedThermalUnit::get_vertical_spacing() {
  if (temp.was_created() == true) {
    return m_Lbz / (m_Mbz - 1.0);
  } else {
    return 0.0;
  }
}


void BedThermalUnit::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  if (temp.was_created()) {
    result.insert(temp.metadata().get_string("short_name"));
  }
}

void BedThermalUnit::define_variables(const std::set<std::string> &vars,
                                                const PIO &nc, IO_Type nctype) {
  if (temp.was_created()) {
    if (set_contains(vars, temp.metadata().get_string("short_name"))) {
      temp.define(nc, nctype);
    }
  }
}

void BedThermalUnit::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (temp.was_created()) {
    if (set_contains(vars, temp.metadata().get_string("short_name"))) {
      temp.write(nc); 
    }
  }
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
void BedThermalUnit::max_timestep(double /*my_t*/, double &my_dt, bool &restrict) {

  if (temp.was_created()) {
    double dzb = this->get_vertical_spacing();
    my_dt = dzb * dzb / (2.0 * bed_D);  // max dt from stability; in seconds
    restrict = true;
  } else {
    my_dt = 0;
    restrict = false;
  }
}


/* FIXME:  the old scheme had better stability properties, as follows:

Because there is no advection, the simplest centered implicit (backward Euler) scheme is easily "bombproof" without choosing \f$\lambda\f$, or other complications.  It has this scaled form,
\anchor bedrockeqn
\f[ -R_b T_{k-1}^{n+1} + \left(1 + 2 R_b\right) T_k^{n+1} - R_b T_{k+1}^{n+1}
         = T_k^n, \tag{bedrockeqn} \f]
where 
  \f[ R_b = \frac{k_b \Delta t}{\rho_b c_b \Delta z^2}. \f]
This is unconditionally stable for a pure bedrock problem, and has a maximum principle, without any further qualification [\ref MortonMayers].

FIXME:  now a trapezoid rule could be used
*/
void BedThermalUnit::update(double my_t, double my_dt) {

  if (temp.was_created() == false) {
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

    if (contiguous == false) {
      throw RuntimeError::formatted("BedThermalUnit::update() requires next update to be contiguous with last;\n"
                                    "  stored:     t = %f s,    dt = %f s\n"
                                    "  desired: my_t = %f s, my_dt = %f s",
                                    m_t,m_dt,my_t,my_dt); }
  }
  // CHECK: is desired time-step too long?
  double my_max_dt;
  bool restrict_dt;
  max_timestep(my_t, my_max_dt, restrict_dt);
  if (restrict_dt && my_max_dt < my_dt) {
     throw RuntimeError("BedThermalUnit::update() thinks you asked for too big a timestep.");
  }

  // o.k., we have checked; we are going to do the desired timestep!
  m_t  = my_t;
  m_dt = my_dt;

  assert(bedtoptemp != NULL);
  assert(ghf != NULL);

  double dzb = this->get_vertical_spacing();
  const int  k0  = m_Mbz - 1;          // Tb[k0] = ice/bed interface temp, at z=0

  const double bed_R  = bed_D * my_dt / (dzb * dzb);

  double *Tbold;
  std::vector<double> Tbnew(m_Mbz);

  IceModelVec::AccessList list;
  list.add(temp);
  list.add(*ghf);
  list.add(*bedtoptemp);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    temp.getInternalColumn(i,j,&Tbold); // Tbold actually points into temp memory
    Tbold[k0] = (*bedtoptemp)(i,j);  // sets Dirichlet explicit-in-time b.c. at top of bedrock column

    const double Tbold_negone = Tbold[1] + 2 * (*ghf)(i,j) * dzb / bed_k;
    Tbnew[0] = Tbold[0] + bed_R * (Tbold_negone - 2 * Tbold[0] + Tbold[1]);
    for (int k = 1; k < k0; k++) { // working upward from base
      Tbnew[k] = Tbold[k] + bed_R * (Tbold[k-1] - 2 * Tbold[k] + Tbold[k+1]);
    }
    Tbnew[k0] = (*bedtoptemp)(i,j);

    temp.setInternalColumn(i,j,&Tbnew[0]); // copy from Tbnew into temp memory
  }
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
void BedThermalUnit::get_upward_geothermal_flux(IceModelVec2S &result) {

  if (!temp.was_created()) {
    result.copy_from(*ghf);
    return;
  }

  double dzb = this->get_vertical_spacing();
  const int  k0  = m_Mbz - 1;  // Tb[k0] = ice/bed interface temp, at z=0

  double *Tb;

  IceModelVec::AccessList list;
  list.add(temp);
  list.add(result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    temp.getInternalColumn(i,j,&Tb);
    if (m_Mbz >= 3) {
      result(i,j) = - bed_k * (3 * Tb[k0] - 4 * Tb[k0-1] + Tb[k0-2]) / (2 * dzb);
    } else {
      result(i,j) = - bed_k * (Tb[k0] - Tb[k0-1]) / dzb;
    }
  }
}

void BedThermalUnit::bootstrap() {

  if (m_Mbz < 2) {
    return;
  }

  verbPrintf(2,m_grid.com,
             "  bootstrapping to fill lithosphere temperatures in bedrock thermal layers,\n"
             "    using provided bedtoptemp and a linear function from provided geothermal flux ...\n");

  double* Tb;
  double dzb = this->get_vertical_spacing();
  const int k0 = m_Mbz-1; // Tb[k0] = ice/bedrock interface temp

  IceModelVec::AccessList list;
  list.add(*bedtoptemp);
  list.add(*ghf);
  list.add(temp);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    temp.getInternalColumn(i,j,&Tb); // Tb points into temp memory
    Tb[k0] = (*bedtoptemp)(i,j);
    for (int k = k0-1; k >= 0; k--) {
      Tb[k] = Tb[k+1] + dzb * (*ghf)(i,j) / bed_k;
    }
  }

  temp.inc_state_counter();     // mark as modified
}


} // end of namespace pism
