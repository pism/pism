// Copyright (C) 2011, 2012 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "LocalInterpCtx.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

bool IceModelVec3BTU::good_init() {
  return ((n_levels >= 2) && (Lbz > 0.0) && (v != PETSC_NULL));
}


PetscErrorCode IceModelVec3BTU::create(IceGrid &mygrid, const char my_short_name[], bool local,
                                      int myMbz, PetscReal myLbz, int stencil_width) {
  PetscErrorCode ierr;
  grid = &mygrid;

  if (!utIsInit()) {
    SETERRQ(grid->com, 1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(grid->com, 2,"IceModelVec3BTU with name='%s' already allocated\n",name.c_str());
  }

  name = my_short_name;

  n_levels = myMbz;
  Lbz = myLbz;
  zlevels.resize(n_levels);
  double dz = Lbz / (myMbz - 1);
  for (int i = 0; i < n_levels; ++i)
    zlevels[i] = -Lbz + i * dz;
  zlevels.back() = 0;

  da_stencil_width = stencil_width;
  ierr = create_2d_da(da, n_levels, da_stencil_width); CHKERRQ(ierr);

  localp = local;
  if (local) {
    ierr = DMCreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  vars[0].init_3d(name, mygrid, zlevels);
  vars[0].dimensions["z"] = "zb";

  map<string,string> &attrs = vars[0].z_attrs;
  attrs["axis"]          = "Z";
  attrs["long_name"]     = "Z-coordinate in bedrock";
  // PROPOSED: attrs["standard_name"] = "projection_z_coordinate_in_lithosphere";
  attrs["units"]         = "m";
  attrs["positive"]      = "up";

  if (!good_init()) {
    SETERRQ1(grid->com, 1,"create() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());  }
  return 0;
}

PetscErrorCode IceModelVec3BTU::get_layer_depth(PetscReal &depth) {
  if (!good_init()) {
    SETERRQ1(grid->com, 1,"get_layer_depth() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());
  }
  depth = Lbz;
  return 0;
}

PetscErrorCode IceModelVec3BTU::get_spacing(PetscReal &dzb) {
  if (!good_init()) {
    SETERRQ1(grid->com, 1,"get_spacing() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());
  }
  dzb = Lbz / (n_levels - 1);
  return 0;
}

PISMBedThermalUnit::PISMBedThermalUnit(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf) {
  bedtoptemp = NULL;
  ghf        = NULL;

  Mbz = 1;
  Lbz = 0;
}

//! \brief Allocate storage for the temperature in the bedrock layer (if there
//! is a bedrock layer).
PetscErrorCode PISMBedThermalUnit::allocate(int my_Mbz, double my_Lbz) {
  PetscErrorCode ierr;

  // to allow multiple calls to init() during the initialization sequence
  if (temp.was_created()) return 0;

  Mbz = my_Mbz;
  Lbz = my_Lbz;

  if ((Lbz <= 0.0) && (Mbz > 1)) {
     SETERRQ(grid.com, 1,"PISMBedThermalUnit can not be created with negative or zero Lbz value\n"
               " and more than one layers\n"); }

  if (Mbz > 1) {
    ierr = temp.create(grid, "litho_temp", false, Mbz, Lbz); CHKERRQ(ierr);
    ierr = temp.set_attrs("model_state",
                          "lithosphere (bedrock) temperature, in PISMBedThermalUnit",
                          "K", ""); CHKERRQ(ierr);
    ierr = temp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  }

  return 0;
}


//! \brief Initialize the bedrock thermal unit.
PetscErrorCode PISMBedThermalUnit::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool i_set, Mbz_set, Lbz_set;
  string input_file;
  grid_info g;

  t = dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2,grid.com,
      "* Initializing the bedrock thermal unit ... setting constants ...\n"); CHKERRQ(ierr);

  // Get pointers to fields owned by IceModel.
  bedtoptemp = dynamic_cast<IceModelVec2S*>(vars.get("bedtoptemp"));
  if (bedtoptemp == NULL) SETERRQ(grid.com, 1, "bedtoptemp is not available");

  ghf = dynamic_cast<IceModelVec2S*>(vars.get("bheatflx"));
  if (ghf == NULL) SETERRQ(grid.com, 2, "bheatflx is not available");

  Mbz = (PetscInt)config.get("grid_Mbz");
  Lbz = (PetscInt)config.get("grid_Lbz");

  // build constant diffusivity for heat equation
  bed_rho = config.get("bedrock_thermal_density");
  bed_c   = config.get("bedrock_thermal_specific_heat_capacity");
  bed_k   = config.get("bedrock_thermal_conductivity");
  bed_D   = bed_k / (bed_rho * bed_c);

  ierr = PetscOptionsBegin(grid.com, "", "PISMBedThermalUnit options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-i", "PISM input file name",
                             input_file, i_set); CHKERRQ(ierr);
    ierr = PISMOptionsInt("-Mbz", "number of levels in bedrock thermal layer", Mbz, Mbz_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-Lbz", "depth (thickness) of bedrock thermal layer", Lbz, Lbz_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (i_set) {
    // If we're initializing from a file we need to get the number of bedrock
    // levels and the depth of the bed thermal layer from it:
    PIO nc(grid.com, grid.rank, "netcdf3");

    ierr = nc.open(input_file, PISM_NOWRITE); CHKERRQ(ierr);
    
    bool exists;
    ierr = nc.inq_var("litho_temp", exists); CHKERRQ(ierr);

    if (exists) {
      ierr = nc.inq_grid_info("litho_temp", g); CHKERRQ(ierr);

      Mbz = g.z_len;
      Lbz = - g.z_min;

      ierr = allocate(Mbz, Lbz); CHKERRQ(ierr);

      int last_record = g.t_len - 1;
      ierr = temp.read(input_file, last_record); CHKERRQ(ierr);
    } else {
      Mbz = 1;
      Lbz = 0;
    }

    ierr = nc.close(); CHKERRQ(ierr);

    ierr = ignore_option(grid.com, "-Mbz"); CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Lbz"); CHKERRQ(ierr);
  } else {
    // Bootstrapping

    if (Mbz_set && Mbz == 1) {
      ierr = ignore_option(grid.com, "-Lbz"); CHKERRQ(ierr);
      Lbz = 0;
    } else if (Mbz_set ^ Lbz_set) {
      PetscPrintf(grid.com, "PISMBedThermalUnit ERROR: please specify both -Mbz and -Lbz.\n");
      PISMEnd();
    }

    ierr = allocate(Mbz, Lbz); CHKERRQ(ierr);

    ierr = bootstrap(); CHKERRQ(ierr);
  }

  // If we're using a minimal model, then we're done:
  if (!temp.was_created()) {
    ierr = verbPrintf(2,grid.com,
      "  minimal model for lithosphere: stored geothermal flux applied to ice base ...\n");
      CHKERRQ(ierr);
      return 0;
  }

  ierr = regrid(); CHKERRQ(ierr);

  return 0;
}


void PISMBedThermalUnit::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  if (temp.was_created()) {
    result[temp.string_attr("short_name")] = temp.get_metadata();
  }
}

PetscErrorCode PISMBedThermalUnit::define_variables(
                         set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  if (temp.was_created()) {
    PetscErrorCode ierr;
    if (set_contains(vars, temp.string_attr("short_name"))) {
      ierr = temp.define(nc, nctype); CHKERRQ(ierr);
    }
  }
  return 0;
}

PetscErrorCode PISMBedThermalUnit::write_variables(set<string> vars, string filename) {
  if (temp.was_created()) {
    PetscErrorCode ierr;
    if (set_contains(vars, temp.string_attr("short_name"))) {
      ierr = temp.write(filename.c_str()); CHKERRQ(ierr); 
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
PetscErrorCode PISMBedThermalUnit::max_timestep(PetscReal /*my_t*/, PetscReal &my_dt, bool &restrict) {

  if (temp.was_created()) {
    PetscReal dzb;
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
PetscErrorCode PISMBedThermalUnit::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if (!temp.was_created())  return 0;  // in this case we are up to date

  // as a derived class of PISMComponent_TS, has t,dt members which keep track
  // of last update time-interval; so we do some checks ...
  // CHECK: has the desired time-interval already been dealt with?
  if ((fabs(my_t - t) < 1e-12) && (fabs(my_dt - dt) < 1e-12))  return 0;

  // CHECK: is the desired time interval a forward step?; backward heat equation not good!
  if (my_dt < 0) {
     SETERRQ(grid.com, 1,"PISMBedThermalUnit::update() does not allow negative timesteps\n"); }
  // CHECK: is desired time-interval equal to [my_t,my_t+my_dt] where my_t = t + dt?
  if ((!gsl_isnan(t)) && (!gsl_isnan(dt))) { // this check should not fire on first use
    bool contiguous = true;

    if (fabs(t + dt) < 1) {
      if ( fabs(my_t - (t + dt)) >= 1e-12 ) // check if the absolute difference is small
        contiguous = false;
    } else {
      if ( fabs(my_t - (t + dt)) / (t + dt) >= 1e-12 ) // check if the relative difference is small
        contiguous = false;
    }

    if (contiguous == false) {
     SETERRQ4(grid.com, 2,"PISMBedThermalUnit::update() requires next update to be contiguous with last;\n"
                "  stored:     t = %f s,    dt = %f s\n"
                "  desired: my_t = %f s, my_dt = %f s\n",
              t,dt,my_t,my_dt); }
  }
  // CHECK: is desired time-step too long?
  PetscScalar my_max_dt;
  bool restrict_dt;
  ierr = max_timestep(my_t, my_max_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt && my_max_dt < my_dt) {
     SETERRQ(grid.com, 3,"PISMBedThermalUnit::update() thinks you asked for too big a timestep\n"); }

  // o.k., we have checked; we are going to do the desired timestep!
  t  = my_t;
  dt = my_dt;

  if (bedtoptemp == NULL)      SETERRQ(grid.com, 5, "bedtoptemp was never initialized");
  if (ghf == NULL)      SETERRQ(grid.com, 6, "bheatflx was never initialized");

  PetscReal dzb;
  temp.get_spacing(dzb);
  const PetscInt  k0  = Mbz - 1;          // Tb[k0] = ice/bed interface temp, at z=0

#if (PISM_DEBUG==1)
  for (PetscInt k = 0; k < Mbz; k++) { // working upward from base
    const PetscReal  z = - Lbz + (double)k * dzb;
    ierr = temp.isLegalLevel(z); CHKERRQ(ierr);
  }
#endif

  const PetscReal bed_R  = bed_D * my_dt / (dzb * dzb);

  PetscScalar *Tbold;
  vector<PetscScalar> Tbnew(Mbz);

  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = ghf->begin_access(); CHKERRQ(ierr);
  ierr = bedtoptemp->begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      ierr = temp.getInternalColumn(i,j,&Tbold); CHKERRQ(ierr); // Tbold actually points into temp memory
      Tbold[k0] = (*bedtoptemp)(i,j);  // sets Dirichlet explicit-in-time b.c. at top of bedrock column

      const PetscReal Tbold_negone = Tbold[1] + 2 * (*ghf)(i,j) * dzb / bed_k;
      Tbnew[0] = Tbold[0] + bed_R * (Tbold_negone - 2 * Tbold[0] + Tbold[1]);
      for (PetscInt k = 1; k < k0; k++) { // working upward from base
        Tbnew[k] = Tbold[k] + bed_R * (Tbold[k-1] - 2 * Tbold[k] + Tbold[k+1]);
      }
      Tbnew[k0] = (*bedtoptemp)(i,j);

      ierr = temp.setInternalColumn(i,j,&Tbnew[0]); CHKERRQ(ierr); // copy from Tbnew into temp memory
    }
  }
  ierr = bedtoptemp->end_access(); CHKERRQ(ierr);
  ierr = ghf->end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);

  return 0;
}


/*! Computes the heat flux from the bedrock thermal layer upward into the
ice/bedrock interface:
  \f[G_0 = -k_b \frac{\partial T_b}{\partial z}\big|_{z=0}.\f]
Uses the second-order finite difference expression
  \f[\frac{\partial T_b}{\partial z}\big|_{z=0} \approx \frac{3 T_b(0) - 4 T_b(-\Delta z) + T_b(-2\Delta z)}{2 \Delta z}\f]
where \f$\Delta z\f$ is the equal spacing in the bedrock.

The above expression only makes sense when \c Mbz = \c temp.n_levels >= 3.
When \c Mbz = 2 we use first-order differencing.  When temp was not created,
the \c Mbz <= 1 cases, we return the stored geothermal flux.
 */
PetscErrorCode PISMBedThermalUnit::get_upward_geothermal_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  if (!temp.was_created()) {
    result.copy_from(*ghf);
    return 0;
  }

  PetscReal dzb;
  temp.get_spacing(dzb);
  const PetscInt  k0  = Mbz - 1;  // Tb[k0] = ice/bed interface temp, at z=0

  PetscScalar *Tb;
  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
      if (Mbz >= 3) {
        result(i,j) = - bed_k * (3 * Tb[k0] - 4 * Tb[k0-1] + Tb[k0-2]) / (2 * dzb);
      } else {
        result(i,j) = - bed_k * (Tb[k0] - Tb[k0-1]) / dzb;
      }
    }
  }
  ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMBedThermalUnit::regrid() {
  PetscErrorCode ierr;
  bool regrid_file_set, regrid_vars_set;
  string regrid_file;
  vector<string> regrid_vars;

  ierr = PetscOptionsBegin(grid.com, "", "PISMBedThermalUnit regridding options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "regridding file name",
                             regrid_file, regrid_file_set); CHKERRQ(ierr);
    ierr = PISMOptionsStringArray("-regrid_vars", "comma-separated list of regridding variables",
                                  "", regrid_vars, regrid_vars_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // stop unless both options are set
  if (! regrid_file_set) return 0;

  // stop if temp was not allocated
  if (! temp.was_created()) return 0;

  set<string> vars;
  for (unsigned int i = 0; i < regrid_vars.size(); ++i)
    vars.insert(regrid_vars[i]);

  // stop if the user did not ask to regrid BTU temperature
  if (regrid_vars_set && ! set_contains(vars, temp.string_attr("short_name")))
    return 0;

  ierr = temp.regrid(regrid_file.c_str(), true); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMBedThermalUnit::bootstrap() {
  PetscErrorCode ierr;

  if (Mbz < 2) return 0;

  ierr = verbPrintf(2,grid.com,
                    "  bootstrapping to fill lithosphere temperatures in bedrock thermal layers,\n"
                    "    using provided bedtoptemp and a linear function from provided geothermal flux ...\n");
  CHKERRQ(ierr);

  PetscScalar* Tb;
  PetscReal dzb;
  temp.get_spacing(dzb);
  const PetscInt k0 = Mbz-1; // Tb[k0] = ice/bedrock interface temp

  ierr = bedtoptemp->begin_access(); CHKERRQ(ierr);
  ierr = ghf->begin_access(); CHKERRQ(ierr);
  ierr = temp.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr); // Tb points into temp memory
      Tb[k0] = (*bedtoptemp)(i,j);
      for (PetscInt k = k0-1; k >= 0; k--) {
        Tb[k] = Tb[k+1] + dzb * (*ghf)(i,j) / bed_k;
      }
    }
  }
  ierr = bedtoptemp->end_access(); CHKERRQ(ierr);
  ierr = ghf->end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);

  return 0;
}

