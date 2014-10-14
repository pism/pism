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

namespace pism {

bool IceModelVec3BTU::good_init() {
  return ((m_n_levels >= 2) && (Lbz > 0.0) && (m_v != NULL));
}


PetscErrorCode IceModelVec3BTU::create(IceGrid &mygrid, const std::string &my_short_name, bool local,
                                      int Mbz, double myLbz, int stencil_width) {
  PetscErrorCode ierr;
  grid = &mygrid;

  if (m_v != NULL) {
    SETERRQ1(grid->com, 2,"IceModelVec3BTU with name='%s' already allocated\n",m_name.c_str());
  }

  m_name = my_short_name;

  m_n_levels = Mbz;
  Lbz = myLbz;
  zlevels.resize(m_n_levels);
  double dz = Lbz / (Mbz - 1);
  for (unsigned int i = 0; i < m_n_levels; ++i)
    zlevels[i] = -Lbz + i * dz;
  zlevels.back() = 0;

  m_da_stencil_width = stencil_width;
  ierr = grid->get_dm(this->m_n_levels, this->m_da_stencil_width, m_da); CHKERRQ(ierr);

  m_has_ghosts = local;
  if (local) {
    ierr = DMCreateLocalVector(*m_da, &m_v); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(*m_da, &m_v); CHKERRQ(ierr);
  }

  m_metadata.resize(m_dof, NCSpatialVariable(grid->get_unit_system()));
  m_metadata[0].init_3d(m_name, mygrid, zlevels);

  NCVariable &z = m_metadata[0].get_z();
  z.set_name("zb");
  z.set_string("axis", "Z");
  z.set_string("long_name", "Z-coordinate in bedrock");
  z.set_string("positive", "up");
  ierr = z.set_units("m"); CHKERRQ(ierr);
  // PROPOSED: attrs["standard_name"] = "projection_z_coordinate_in_lithosphere";

  if (!good_init()) {
    SETERRQ1(grid->com, 1,"create() says IceModelVec3BTU with name %s was not properly created\n",
             m_name.c_str());  }
  return 0;
}

PetscErrorCode IceModelVec3BTU::get_layer_depth(double &depth) {
  if (!good_init()) {
    SETERRQ1(grid->com, 1,"get_layer_depth() says IceModelVec3BTU with name %s was not properly created\n",
             m_name.c_str());
  }
  depth = Lbz;
  return 0;
}

PetscErrorCode IceModelVec3BTU::get_spacing(double &dzb) {
  if (!good_init()) {
    SETERRQ1(grid->com, 1,"get_spacing() says IceModelVec3BTU with name %s was not properly created\n",
             m_name.c_str());
  }
  dzb = Lbz / (m_n_levels - 1);
  return 0;
}

BedThermalUnit::BedThermalUnit(IceGrid &g, const Config &conf)
    : Component_TS(g, conf) {
  bedtoptemp = NULL;
  ghf        = NULL;

  // build constant diffusivity for heat equation
  bed_rho = config.get("bedrock_thermal_density");
  bed_c   = config.get("bedrock_thermal_specific_heat_capacity");
  bed_k   = config.get("bedrock_thermal_conductivity");
  bed_D   = bed_k / (bed_rho * bed_c);

  Mbz = (int)config.get("grid_Mbz");
  Lbz = (int)config.get("grid_Lbz");
  m_input_file.clear();

  PetscErrorCode ierr = allocate(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

//! \brief Allocate storage for the temperature in the bedrock layer (if there
//! is a bedrock layer).
PetscErrorCode BedThermalUnit::allocate() {
  PetscErrorCode ierr;
  bool i_set, Mbz_set, Lbz_set;
  grid_info g;

  ierr = PetscOptionsBegin(grid.com, "", "BedThermalUnit options", ""); CHKERRQ(ierr);
  {
    ierr = OptionsString("-i", "PISM input file name",
                             m_input_file, i_set); CHKERRQ(ierr);

    int tmp = Mbz;
    ierr = OptionsInt("-Mbz", "number of levels in bedrock thermal layer",
                          tmp, Mbz_set); CHKERRQ(ierr);
    Mbz = tmp;

    ierr = OptionsReal("-Lbz", "depth (thickness) of bedrock thermal layer, in meters",
                           Lbz, Lbz_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);


  if (i_set) {
    ierr = ignore_option(grid.com, "-Mbz"); CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Lbz"); CHKERRQ(ierr);

    // If we're initializing from a file we need to get the number of bedrock
    // levels and the depth of the bed thermal layer from it:
    PIO nc(grid, "guess_mode");

    ierr = nc.open(m_input_file, PISM_READONLY); CHKERRQ(ierr);

    bool exists;
    ierr = nc.inq_var("litho_temp", exists); CHKERRQ(ierr);

    if (exists) {
      ierr = nc.inq_grid_info("litho_temp", grid.periodicity, g); CHKERRQ(ierr);

      Mbz = g.z_len;
      Lbz = -g.z_min;
    } else {
      // override values we got using config.get() in the constructor
      Mbz = 1;
      Lbz = 0;
    }

    ierr = nc.close(); CHKERRQ(ierr);
  } else {
    // Bootstrapping

    if (Mbz_set && Mbz == 1) {
      ierr = ignore_option(grid.com, "-Lbz"); CHKERRQ(ierr);
      Lbz = 0;
    } else if (Mbz_set ^ Lbz_set) {
      PetscPrintf(grid.com, "BedThermalUnit ERROR: please specify both -Mbz and -Lbz.\n");
      PISMEnd();
    }
  }

  // actual allocation:
  if ((Lbz <= 0.0) && (Mbz > 1)) {
    SETERRQ(grid.com, 1,"BedThermalUnit can not be created with negative or zero Lbz value\n"
            " and more than one layers\n"); }

  if (Mbz > 1) {
    ierr = temp.create(grid, "litho_temp", false, Mbz, Lbz); CHKERRQ(ierr);
    ierr = temp.set_attrs("model_state",
                          "lithosphere (bedrock) temperature, in BedThermalUnit",
                          "K", ""); CHKERRQ(ierr);
    temp.metadata().set_double("valid_min", 0.0);
  }

  return 0;
}


//! \brief Initialize the bedrock thermal unit.
PetscErrorCode BedThermalUnit::init(Vars &vars, bool &bootstrapping_needed) {
  PetscErrorCode ierr;
  grid_info g;

  // first assume that we don't need to bootstrap
  bootstrapping_needed = false;

  // store the current "revision number" of the temperature field
  int temp_revision = temp.get_state_counter();

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2,grid.com,
      "* Initializing the bedrock thermal unit... setting constants...\n"); CHKERRQ(ierr);

  // Get pointers to fields owned by IceModel.
  bedtoptemp = vars.get_2d_scalar("bedtoptemp");
  ghf = vars.get_2d_scalar("bheatflx");

  // If we're using a minimal model, then we're done:
  if (!temp.was_created()) {
    ierr = verbPrintf(2,grid.com,
      "  minimal model for lithosphere: stored geothermal flux applied to ice base ...\n");
      CHKERRQ(ierr);
      return 0;
  }

  if (m_input_file.empty() == false) {
    PIO nc(grid, "guess_mode");
    bool exists;
    unsigned int n_records;

    ierr = nc.open(m_input_file, PISM_READONLY); CHKERRQ(ierr);
    ierr = nc.inq_var("litho_temp", exists); CHKERRQ(ierr);

    if (exists) {
      ierr = nc.inq_nrecords("litho_temp", "", n_records); CHKERRQ(ierr);

      const unsigned int last_record = n_records - 1;
      ierr = temp.read(m_input_file, last_record); CHKERRQ(ierr);
    }

    ierr = nc.close(); CHKERRQ(ierr);
  }

  if (temp.was_created() == true) {
    ierr = regrid("BedThermalUnit", &temp, REGRID_WITHOUT_REGRID_VARS); CHKERRQ(ierr);
  }

  if (temp.get_state_counter() == temp_revision) {
    bootstrapping_needed = true;
  }

  return 0;
}

/** Returns the vertical spacing used by the bedrock grid.
 *
 * Special case: returns 0 if the bedrock thermal layer has thickness
 * zero.
 */
double BedThermalUnit::get_vertical_spacing() {
  if (temp.was_created() == true) {
    double dzb = 0.0;
    temp.get_spacing(dzb);
    return dzb;
  } else {
    return 0.0;
  }
}


void BedThermalUnit::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  if (temp.was_created()) {
    result.insert(temp.metadata().get_string("short_name"));
  }
}

PetscErrorCode BedThermalUnit::define_variables(const std::set<std::string> &vars,
                                                const PIO &nc, IO_Type nctype) {
  if (temp.was_created()) {
    PetscErrorCode ierr;
    if (set_contains(vars, temp.metadata().get_string("short_name"))) {
      ierr = temp.define(nc, nctype); CHKERRQ(ierr);
    }
  }
  return 0;
}

PetscErrorCode BedThermalUnit::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (temp.was_created()) {
    PetscErrorCode ierr;
    if (set_contains(vars, temp.metadata().get_string("short_name"))) {
      ierr = temp.write(nc); CHKERRQ(ierr); 
    }
  }
  return 0;
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
PetscErrorCode BedThermalUnit::max_timestep(double /*my_t*/, double &my_dt, bool &restrict) {

  if (temp.was_created()) {
    double dzb;
    temp.get_spacing(dzb);
    my_dt = dzb * dzb / (2.0 * bed_D);  // max dt from stability; in seconds
    restrict = true;
  } else {
    my_dt = 0;
    restrict = false;
  }
  return 0;
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
PetscErrorCode BedThermalUnit::update(double my_t, double my_dt) {
  PetscErrorCode ierr;

  if (temp.was_created() == false)
    return 0;  // in this case we are up to date

  // as a derived class of Component_TS, has t,dt members which keep track
  // of last update time-interval; so we do some checks ...
  // CHECK: has the desired time-interval already been dealt with?
  if ((fabs(my_t - m_t) < 1e-12) && (fabs(my_dt - m_dt) < 1e-12))
    return 0;

  // CHECK: is the desired time interval a forward step?; backward heat equation not good!
  if (my_dt < 0) {
     SETERRQ(grid.com, 1,"BedThermalUnit::update() does not allow negative timesteps\n"); }
  // CHECK: is desired time-interval equal to [my_t,my_t+my_dt] where my_t = t + dt?
  if ((!gsl_isnan(m_t)) && (!gsl_isnan(m_dt))) { // this check should not fire on first use
    bool contiguous = true;

    if (fabs(m_t + m_dt) < 1) {
      if (fabs(my_t - (m_t + m_dt)) >= 1e-12) // check if the absolute difference is small
        contiguous = false;
    } else {
      if (fabs(my_t - (m_t + m_dt)) / (m_t + m_dt) >= 1e-12) // check if the relative difference is small
        contiguous = false;
    }

    if (contiguous == false) {
     SETERRQ4(grid.com, 2,"BedThermalUnit::update() requires next update to be contiguous with last;\n"
                "  stored:     t = %f s,    dt = %f s\n"
                "  desired: my_t = %f s, my_dt = %f s\n",
              m_t,m_dt,my_t,my_dt); }
  }
  // CHECK: is desired time-step too long?
  double my_max_dt;
  bool restrict_dt;
  ierr = max_timestep(my_t, my_max_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt && my_max_dt < my_dt) {
     SETERRQ(grid.com, 3,"BedThermalUnit::update() thinks you asked for too big a timestep\n"); }

  // o.k., we have checked; we are going to do the desired timestep!
  m_t  = my_t;
  m_dt = my_dt;

  assert(bedtoptemp != NULL);
  assert(ghf != NULL);

  double dzb;
  temp.get_spacing(dzb);
  const int  k0  = Mbz - 1;          // Tb[k0] = ice/bed interface temp, at z=0

#if (PISM_DEBUG==1)
  for (unsigned int k = 0; k < Mbz; k++) { // working upward from base
    const double  z = - Lbz + (double)k * dzb;
    ierr = temp.isLegalLevel(z); CHKERRQ(ierr);
  }
#endif

  const double bed_R  = bed_D * my_dt / (dzb * dzb);

  double *Tbold;
  std::vector<double> Tbnew(Mbz);

  IceModelVec::AccessList list;
  list.add(temp);
  list.add(*ghf);
  list.add(*bedtoptemp);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = temp.getInternalColumn(i,j,&Tbold); CHKERRQ(ierr); // Tbold actually points into temp memory
    Tbold[k0] = (*bedtoptemp)(i,j);  // sets Dirichlet explicit-in-time b.c. at top of bedrock column

    const double Tbold_negone = Tbold[1] + 2 * (*ghf)(i,j) * dzb / bed_k;
    Tbnew[0] = Tbold[0] + bed_R * (Tbold_negone - 2 * Tbold[0] + Tbold[1]);
    for (int k = 1; k < k0; k++) { // working upward from base
      Tbnew[k] = Tbold[k] + bed_R * (Tbold[k-1] - 2 * Tbold[k] + Tbold[k+1]);
    }
    Tbnew[k0] = (*bedtoptemp)(i,j);

    ierr = temp.setInternalColumn(i,j,&Tbnew[0]); CHKERRQ(ierr); // copy from Tbnew into temp memory
  }

  return 0;
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
PetscErrorCode BedThermalUnit::get_upward_geothermal_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (!temp.was_created()) {
    result.copy_from(*ghf);
    return 0;
  }

  double dzb;
  temp.get_spacing(dzb);
  const int  k0  = Mbz - 1;  // Tb[k0] = ice/bed interface temp, at z=0

  double *Tb;

  IceModelVec::AccessList list;
  list.add(temp);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
    if (Mbz >= 3) {
      result(i,j) = - bed_k * (3 * Tb[k0] - 4 * Tb[k0-1] + Tb[k0-2]) / (2 * dzb);
    } else {
      result(i,j) = - bed_k * (Tb[k0] - Tb[k0-1]) / dzb;
    }
  }

  return 0;
}

PetscErrorCode BedThermalUnit::bootstrap() {
  PetscErrorCode ierr;

  if (Mbz < 2) return 0;

  ierr = verbPrintf(2,grid.com,
                    "  bootstrapping to fill lithosphere temperatures in bedrock thermal layers,\n"
                    "    using provided bedtoptemp and a linear function from provided geothermal flux ...\n");
  CHKERRQ(ierr);

  double* Tb;
  double dzb;
  temp.get_spacing(dzb);
  const int k0 = Mbz-1; // Tb[k0] = ice/bedrock interface temp

  IceModelVec::AccessList list;
  list.add(*bedtoptemp);
  list.add(*ghf);
  list.add(temp);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr); // Tb points into temp memory
    Tb[k0] = (*bedtoptemp)(i,j);
    for (int k = k0-1; k >= 0; k--) {
      Tb[k] = Tb[k+1] + dzb * (*ghf)(i,j) / bed_k;
    }
  }

  temp.inc_state_counter();     // mark as modified

  return 0;
}


} // end of namespace pism
