// Copyright (C) 2011 Ed Bueler
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


bool IceModelVec3BTU::good_init() {
  return ((n_levels >= 3) && (Lbz > 0.0) && (v != PETSC_NULL));
}


PetscErrorCode  IceModelVec3BTU::create(IceGrid &mygrid, const char my_short_name[], bool local,
                                        int myMbz, PetscReal myLbz, int stencil_width) {
  PetscErrorCode ierr;
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec3BTU with name='%s' already allocated\n",name.c_str());
  }

  grid = &mygrid;
  name = my_short_name;

  dims = GRID_3D_BEDROCK;  // FIXME:  not sure what this means ...
  n_levels = myMbz;
  Lbz = myLbz;

  da_stencil_width = stencil_width;
  ierr = create_2d_da(da, n_levels, da_stencil_width); CHKERRQ(ierr);

  localp = local;
  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  vars[0].init(name, mygrid, dims);

  if (!good_init()) {
    SETERRQ1(1,"create() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());
  }
  return 0;
}


PetscErrorCode IceModelVec3BTU::get_levels(PetscInt &levels) {
  if (!good_init()) {
    SETERRQ1(1,"get_levels() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());
  }
  levels = n_levels;
  return 0;
}


PetscErrorCode IceModelVec3BTU::get_layer_depth(PetscReal &depth) {
  if (!good_init()) {
    SETERRQ1(1,"get_layer_depth() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());
  }
  depth = Lbz;
  return 0;
}


PetscErrorCode IceModelVec3BTU::get_spacing(PetscReal &dzb) {
  if (!good_init()) {
    SETERRQ1(1,"get_spacing() says IceModelVec3BTU with name %s was not properly created\n",
             name.c_str());
  }
  dzb = Lbz / (n_levels - 1);
  return 0;
}


PetscErrorCode IceModelVec3BTU::stopIfNotLegalLevel(PetscScalar z) {
  if (z < -Lbz) {
    SETERRQ3(1,
       "level z = %10.8f is below bottom of bedrock thermal layer at -Lbz = %10.8f;\n"
       "IceModelVec3BTU has name='%s'; ENDING!\n",
       z,-Lbz,name.c_str());
  }
  if (z > 0.0) {
    SETERRQ2(2,"level z = %10.8f is above top of bedrock at z=0;\n"
               " IceModelVec3BTU has name='%s'; ENDING!\n",
               z,name.c_str());
  }
  return 0;
}


PISMBedThermalUnit::PISMBedThermalUnit(IceGrid &g, EnthalpyConverter &e, 
                                       const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf), EC(e) {

  enthalpy = NULL;
  thk      = NULL;
  ghf      = NULL;
  if (allocate()) {
    verbPrintf(1,g.com, "PISMBedThermalUnit::allocate() returned nonzero\n");
    PISMEnd();
  }
}


PetscErrorCode PISMBedThermalUnit::allocate() {
  PetscErrorCode ierr;

  // FIXME: this "raw" way of getting the size should be reevaluated
  bool flag;
  PetscInt Mbz;
  ierr = PISMOptionsInt("-Mbz", "number of levels in bedrock thermal layer", Mbz, flag); CHKERRQ(ierr);
  if (!flag) {
     SETERRQ(1,"option -Mbz was not set so PISMBedThermalUnit cannot be allocated\n"); }
  if (Mbz < 3) {
     SETERRQ(2,"Mbz < 3, but PISMBedThermalUnit creation requires at least 3 levels\n"); }

  PetscReal Lbz;
  ierr = PISMOptionsReal("-Lbz", "depth (thickness) of bedrock thermal layer", Lbz, flag); CHKERRQ(ierr);
  if (!flag) {
     SETERRQ(3,"option -Lbz was not set so PISMBedThermalUnit cannot be allocated\n"); }
  if (Lbz <= 0.0) {
     SETERRQ(4,"PISMBedThermalUnit can not be created with negative or zero Lbz value\n"); }

  ierr = temp.create(grid, "litho_temp", false, Mbz, Lbz); CHKERRQ(ierr);
  ierr = temp.set_attrs("model_state",
                        "lithosphere (bedrock) temperature, in PISMBedThermalUnit",
		        "K", ""); CHKERRQ(ierr);
  ierr = temp.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = ice_base_temp.create(grid, "btu_ice_base_temp", false); CHKERRQ(ierr);
  ierr = ice_base_temp.set_attrs("internal",
    "temperature at base of ice for duration of timestep, in PISMBedThermalUnit",
    "K", ""); CHKERRQ(ierr);
  ierr = ice_base_temp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  
  return 0;
}


/*! this init method is really a bootstrap method.

Note \c variables.get("enthalpy") and  \c variables.get("thk") and 
\c variables.get("bheatflx") must return valid pointers.
 */
PetscErrorCode PISMBedThermalUnit::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com,
      "* Initializing the bedrock thermal unit ... setting constants ...\n"); CHKERRQ(ierr);

  // build constant diffusivity for heat equation
  bed_rho = config.get("bedrock_thermal_density");
  bed_c   = config.get("bedrock_thermal_specific_heat_capacity");
  bed_k   = config.get("bedrock_thermal_conductivity");
  bed_D   = bed_k / (bed_rho * bed_c);

  ierr = verbPrintf(2,grid.com,
      "  bootstrapping using z=0 enthalpy values and linear function from geothermal flux ...\n");
      CHKERRQ(ierr);

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(2, "enthalpy is not available");
  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(3, "thk is not available");
  ghf = dynamic_cast<IceModelVec2S*>(vars.get("bheatflx"));
  if (ghf == NULL) SETERRQ(4, "btu_bheatflx is not available");

  // first job: fill ice_base_temp; FIXME: is correct in floating case? ice-free case?; FIXME we need mask
  ierr = enthalpy->getHorSlice(ice_base_temp, 0.0); CHKERRQ(ierr);

  PetscScalar* Tb;
  PetscReal dzb;
  PetscInt  Mbz;
  temp.get_levels(Mbz);
  temp.get_spacing(dzb);
  const PetscInt k0 = Mbz-1; // Tb[k0] = ice/bedrock interface temp

  ierr = thk->begin_access(); CHKERRQ(ierr);
  ierr = ghf->begin_access(); CHKERRQ(ierr);
  ierr = ice_base_temp.begin_access(); CHKERRQ(ierr);
  ierr = temp.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      const PetscReal pressure = EC.getPressureFromDepth((*thk)(i,j));
      PetscReal base_temp;
      ierr = EC.getAbsTemp(ice_base_temp(i,j), pressure, base_temp); CHKERRQ(ierr);
      ice_base_temp(i,j) = base_temp;
      ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr); // Tb points into temp memory
      Tb[k0] = ice_base_temp(i,j);
      for (PetscInt k = k0-1; k >= 0; k--) {
        Tb[k] = Tb[k+1] + dzb * (*ghf)(i,j) / bed_k;
      }
    }
  }
  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = ghf->end_access(); CHKERRQ(ierr);
  ierr = ice_base_temp.end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);

  return 0;
}


// FIXME:  init from file code is needed when running with IceModel, at least


void PISMBedThermalUnit::add_vars_to_output(string keyword, set<string> &result) {
  result.insert(temp.string_attr("short_name"));
  if (keyword == "big") {
    result.insert(ice_base_temp.string_attr("short_name"));
  }
}


PetscErrorCode PISMBedThermalUnit::define_variables(
                         set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, temp.string_attr("short_name"))) {
    ierr = temp.define(nc, nctype); CHKERRQ(ierr); 
  }  
  if (set_contains(vars, ice_base_temp.string_attr("short_name"))) {
    ierr = ice_base_temp.define(nc, nctype); CHKERRQ(ierr); 
  }  
  return 0;
}


PetscErrorCode PISMBedThermalUnit::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  if (set_contains(vars, temp.string_attr("short_name"))) {
    ierr = temp.write(filename.c_str()); CHKERRQ(ierr); 
  }  
  if (set_contains(vars, ice_base_temp.string_attr("short_name"))) {
    ierr = ice_base_temp.write(filename.c_str()); CHKERRQ(ierr); 
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
 */
PetscErrorCode PISMBedThermalUnit::max_timestep(PetscReal /*t_years*/, PetscReal &dt_years) {

  PetscReal dzb;
  temp.get_spacing(dzb);
  dt_years = dzb * dzb / (2.0 * bed_D);  // max dt from stability; in seconds
  dt_years /= secpera;
  return 0;
}


PetscErrorCode PISMBedThermalUnit::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  // as a derived class of PISMComponent_TS, has t,dt members which keep track
  // of last update time-interval; so we do some checks ...
  // CHECK: has the desired time-interval already been dealt with?
  if ((fabs(t_years - t) < 1e-12) && (fabs(dt_years - dt) < 1e-12))
    return 0;
  // CHECK: is the desired time interval a forward step?; backward heat equation not good!
  if (dt_years < 0) {
     SETERRQ(1,"PISMBedThermalUnit::update() does not allow negative timesteps\n"); }
  // CHECK: is desired time-interval equal to [t_years,t_years+dt_years] where t_years = t + dt?
  if ((!gsl_isnan(t)) && (!gsl_isnan(dt))) { // this check should not fire on first use
    if (!(fabs(t_years - (t + dt) < 1e-12))) {
     SETERRQ3(2,"PISMBedThermalUnit::update() requires next update to be contiguous with last;\n"
                "situation: t = %f a, dt = %f a, t_years = %f a\n",
              t,dt,t_years); }
  }
  // CHECK: is desired time-step too long?
  PetscScalar mydtyears;
  ierr = max_timestep(t_years,mydtyears); CHKERRQ(ierr);
  if (mydtyears < dt_years) {
     SETERRQ(3,"PISMBedThermalUnit::update() thinks you asked for too big a timestep\n"); }

  // o.k., we have checked; we are going to do the desired timestep!
  t  = t_years;
  dt = dt_years;

  if (enthalpy == NULL) SETERRQ(4, "enthalpy was never initialized");
  if (thk == NULL)      SETERRQ(5, "thk was never initialized");
  if (ghf == NULL)      SETERRQ(6, "bheatflx was never initialized");

  // first job: fill ice_base_temp; FIXME: is correct in floating case? ice-free case?; FIXME we need mask
  ierr = enthalpy->getHorSlice(ice_base_temp, 0.0); CHKERRQ(ierr);
  ierr = ice_base_temp.begin_access(); CHKERRQ(ierr);
  ierr = thk->begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      const PetscReal pressure = EC.getPressureFromDepth((*thk)(i,j));
      ierr = EC.getAbsTemp(ice_base_temp(i,j), pressure, ice_base_temp(i,j));
               CHKERRQ(ierr);
    }
  }
  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = ice_base_temp.end_access(); CHKERRQ(ierr);
  // no need to communicate just-filled ice_base_temp

  PetscReal dzb;
  PetscInt  Mbz;
  temp.get_levels(Mbz);
  temp.get_spacing(dzb);
  const PetscInt  k0  = Mbz - 1;          // Tb[k0] = ice/bed interface temp, at z=0

#ifdef PISM_DEBUG
  PetscReal Lbz;
  temp.get_layer_depth(Lbz);
  for (PetscInt k = 0; k < Mbz; k++) { // working upward from base
    const PetscReal  z = - Lbz + (double)k * dzb;
    ierr = temp.stopIfNotLegalLevel(z); CHKERRQ(ierr);
  }
#endif

  const PetscReal bed_R  = bed_D * (dt_years * secpera) / (dzb * dzb);

  PetscScalar *Tbold, *Tbnew;
  Tbnew = new PetscScalar[Mbz];

  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = ghf->begin_access(); CHKERRQ(ierr);
  ierr = ice_base_temp.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      ierr = temp.getInternalColumn(i,j,&Tbold); CHKERRQ(ierr); // Tbold actually points into temp memory
      Tbold[k0] = ice_base_temp(i,j);  // sets Dirichlet explicit-in-time b.c. at top of bedrock column

      const PetscReal Tbold_negone = Tbold[1] + 2 * (*ghf)(i,j) * dzb / bed_k;
      Tbnew[0] = Tbold[0] + bed_R * (Tbold_negone - 2 * Tbold[0] + Tbold[1]);
      for (PetscInt k = 1; k < k0; k++) { // working upward from base
        Tbnew[k] = Tbold[k] + bed_R * (Tbold[k-1] - 2 * Tbold[k] + Tbold[k+1]);
      }
      Tbnew[k0] = ice_base_temp(i,j);

      ierr = temp.setInternalColumn(i,j,Tbnew); CHKERRQ(ierr); // copy from Tbnew into temp memory
    }
  }
  ierr = ice_base_temp.end_access(); CHKERRQ(ierr);
  ierr = ghf->end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);

  delete [] Tbnew;

  return 0;
}


/*! Computes the heat flux from the bedrock thermal layer upward into the
ice/bedrock interface:
  \f[G_0 = -k_b \frac{\partial T_b}{\partial z}\big|_{z=0}.\f]
Uses the second-order finite difference expression
  \f[\frac{\partial T_b}{\partial z}\big|_{z=0} \approx \frac{3 T_b(0) - 4 T_b(-\Delta z) + T_b(-2\Delta z)}{2 \Delta z}\f]
where \f$\Delta z\f$ is the equal spacing in the bedrock.  This expression
explains why \c n_levels >= 3 is required.
 */
PetscErrorCode PISMBedThermalUnit::get_upward_geothermal_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  PetscReal dzb;
  PetscInt  Mbz;
  temp.get_levels(Mbz);
  temp.get_spacing(dzb);
  const PetscInt  k0  = Mbz - 1;  // Tb[k0] = ice/bed interface temp, at z=0

  PetscScalar *Tb;
  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
      result(i,j) = - bed_k * (3 * Tb[k0] - 4 * Tb[k0-1] + Tb[k0-2]) / (2 * dzb);  
    }
  }
  ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

// FIXME: a derived class version of PISMBedThermalUnit::get_upward_geothermal_flux()
//        will be simpler: it will just report the values of the stored geothermal flux;
// FIXME: for such a derived class, update() will be null and always succeed

