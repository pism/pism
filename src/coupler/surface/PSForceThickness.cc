// Copyright (C) 2011, 2012 PISM Authors
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

#include "PSForceThickness.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "PIO.hh"

///// "Force-to-thickness" mechanism

void PSForceThickness::attach_atmosphere_model(PISMAtmosphereModel *input) {
  input_model->attach_atmosphere_model(input);
}

PetscErrorCode PSForceThickness::init(PISMVars &vars) {
  PetscErrorCode ierr;
  char fttfile[PETSC_MAX_PATH_LEN] = "";
  PetscBool opt_set;
  PetscScalar fttalpha;
  PetscBool  fttalphaSet;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Surface model forcing", ""); CHKERRQ(ierr);

  ierr = PetscOptionsString("-force_to_thk",
			    "Specifies the target thickness file for the force-to-thickness mechanism",
			    "", "",
			    fttfile, PETSC_MAX_PATH_LEN, &opt_set); CHKERRQ(ierr);

  if (!opt_set) {
    ierr = PetscPrintf(grid.com,
      "ERROR: surface model forcing requires the -force_to_thk option.\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = PetscOptionsReal("-force_to_thk_alpha",
			  "Specifies the force-to-thickness alpha value in per-year units",
			  "", convert(alpha,"yr-1","s-1"),
			  &fttalpha, &fttalphaSet); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "* Initializing force-to-thickness mass-balance modifier...\n"); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) SETERRQ(grid.com, 1, "ERROR: land_ice_thickness is not available");

  ierr = target_thickness.create(grid, "thk", false); CHKERRQ(ierr); // name to read by
  ierr = target_thickness.set_attrs( // set attributes for the read stage; see below for reset
     "diagnostic",
     "target thickness for force-to-thickness mechanism (hit this at end of run)",
     "m",
     "land_ice_thickness"); CHKERRQ(ierr); // standard_name *to read by*
  target_thickness.write_in_glaciological_units = true;

  ierr = ftt_mask.create(grid, "ftt_mask", false); CHKERRQ(ierr);
  ierr = ftt_mask.set_attrs(
     "diagnostic",
     "mask specifying where to apply the force-to-thickness mechanism",
     "", ""); CHKERRQ(ierr); // no units and no standard name
  ierr = ftt_mask.set(1.0); CHKERRQ(ierr); // default: applied in whole domain
  ftt_mask.write_in_glaciological_units = true;

  input_file = fttfile;

  // determine exponential rate alpha from user option or from factor; option
  // is given in a^{-1}
  if (fttalphaSet == PETSC_TRUE) {
    ierr = verbPrintf(3, grid.com, "    option -force_to_thk_alpha seen\n");
       CHKERRQ(ierr);
    alpha = convert(fttalpha,"yr-1","s-1");
  }

  ierr = verbPrintf(2, grid.com,
		    "    alpha = %.6f a-1 for -force_to_thk mechanism\n",
		    convert(alpha,"s-1","yr-1")); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // fttfile now contains name of -force_to_thk file; now check
  // it is really there; and regrid the target thickness
  PIO nc(grid.com, grid.rank, "netcdf3");
  bool mask_exists = false;
  ierr = nc.open(fttfile, PISM_NOWRITE); CHKERRQ(ierr);
  ierr = nc.inq_var("ftt_mask", mask_exists); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "    reading target thickness 'thk' from %s ...\n"
		    "    (this field will appear in output file as 'ftt_target_thk')\n",
		    fttfile); CHKERRQ(ierr);
  ierr = target_thickness.regrid(fttfile, true); CHKERRQ(ierr);

  if (mask_exists) {
    ierr = verbPrintf(2, grid.com,
                      "    reading force-to-thickness mask 'ftt_mask' from %s ...\n", fttfile); CHKERRQ(ierr);
    ierr = ftt_mask.regrid(fttfile, true); CHKERRQ(ierr);
  }

  // reset name to avoid confusion; set attributes again to overwrite "read by" choices above
  ierr = target_thickness.set_name("ftt_target_thk"); CHKERRQ(ierr);
  ierr = target_thickness.set_attrs(
    "diagnostic",
    "target thickness for force-to-thickness mechanism (wants to hit this at end of run)",
    "m",
    "");  // no CF standard_name, to put it mildly
    CHKERRQ(ierr);
  target_thickness.write_in_glaciological_units = true;

  climatic_mass_balance.init_2d("climatic_mass_balance", grid);
  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                  "ice-equivalent surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                  "land_ice_surface_specific_mass_balance");
  ierr = climatic_mass_balance.set_units("m s-1"); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}

/*!
If \c -force_to_thk \c foo.nc is in use then vthktarget will have a target ice thickness
map.  Let \f$H_{\text{tar}}\f$ be this target thickness,
and let \f$H\f$ be the current model thickness.  Recall that the mass continuity
equation solved by IceModel::massContExplicitStep() is
  \f[ \frac{\partial H}{\partial t} = M - S - \nabla\cdot \mathbf{q} \f]
and that this procedure is supposed to produce \f$M\f$.
In this context, the semantics of \c -force_to_thk are that \f$M\f$ is modified
by a multiple of the difference between the target thickness and the current thickness.
In particular, the \f$\Delta M\f$ that is produced here is
  \f[\Delta M = \alpha (H_{\text{tar}} - H)\f]
where \f$\alpha>0\f$ is determined below.  Note \f$\Delta M\f$ is positive in
areas where \f$H_{\text{tar}} > H\f$, so we are adding mass there, and we are ablating
in the other case.

Let \f$t_s\f$ be the start time and \f$t_e\f$ the end time for the run.
Without flow or basal mass balance, or any surface mass balance other than the
\f$\Delta M\f$ computed here, we are solving
  \f[ \frac{\partial H}{\partial t} = \alpha (H_{\text{tar}} - H) \f]
Let's assume \f$H(t_s)=H_0\f$.  This initial value problem has solution
\f$H(t) = H_{\text{tar}} + (H_0 - H_{\text{tar}}) e^{-\alpha (t-t_s)}\f$
and so
  \f[ H(t_e) = H_{\text{tar}} + (H_0 - H_{\text{tar}}) e^{-\alpha (t_e-t_s)} \f]

The constant \f$\alpha\f$ has a default value \c pism_config:force_to_thickness_alpha.

The next example uses files generated from the EISMINT-Greenland experiment;
see the corresponding chapter of the User's Manual.

Suppose we regard the SSL2 run as a spin-up to reach a better temperature field.
It is a spinup in which the surface was allowed to evolve.  Assume the
early file \c green20km_y1.nc has the target thickness, because it essentially
has the input thickness.  This script adds a 500 a run, to finalize the spinup,
in which the ice sheet geometry goes from the the thickness values in
\c green_ssl2_110ka.nc to values very close to those in \c green20km_y1.nc:
\code
#!/bin/bash

NN=8  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "test_ftt.sh 8" then NN = 8
  NN="$1"
fi

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var DO is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi

# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=""    # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pismr -cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi


PISM="${PISM_PREFIX}${PISM_EXEC}"

cmd="$PISM_MPIDO $NN $PISM -ys -1000.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd -pdd_fausto \
    -o no_force.nc -ts_file ts_no_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -1000.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc \
  -o default_force.nc -ts_file ts_default_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -1000.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc -force_to_thk_alpha 0.005 \
    -o weak_force.nc -ts_file ts_weak_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd


cmd="$PISM_MPIDO $NN $PISM -ys -1000.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc -force_to_thk_alpha 0.05 \
    -o strong_force.nc -ts_file ts_strong_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd

\endcode
The script also has a run with no forcing, one with forcing at a lower alpha value,
a factor of five smaller than the default, and one with a forcing at a higher alpha value, a factor of five higher.

As shown below, the time series for \c ivol and \c maximum_diffusivity in the
above time series files show that the force-to-thickness mechanism is forcing
a system with negative feedback.

\image html ivol_force_to_thk.png "\b Volume results from the -force_to_thk mechanism."
\anchor ivol_force_to_thk

\image html diffusivity_force_to_thk.png "\b Maximum diffusivity results from the -force_to_thk mechanism."
\anchor diffusivity_force_to_thk

 */
PetscErrorCode PSForceThickness::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  // get the surface mass balance result from the next level up
  ierr = input_model->ice_surface_mass_flux(result); CHKERRQ(ierr);

  ierr = verbPrintf(5, grid.com,
     "    updating surface mass balance using -force_to_thk mechanism ...");
     CHKERRQ(ierr);

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = target_thickness.begin_access(); CHKERRQ(ierr);
  ierr = ftt_mask.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (ftt_mask(i,j) > 0.5) {
        result(i,j) += alpha * (target_thickness(i,j) - H[i][j]);
      }
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = target_thickness.end_access(); CHKERRQ(ierr);
  ierr = ftt_mask.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  // no communication needed

  return 0;
}

//! Does not modify ice surface temperature.
PetscErrorCode PSForceThickness::ice_surface_temperature(IceModelVec2S &result) {
  return input_model->ice_surface_temperature(result);
}

/*!
The timestep restriction is, by direct analogy, the same as for
   \f[\frac{dy}{dt} = - \alpha y\f]
with explicit (forward) Euler.  If \f$\Delta t\f$ is the time step then Euler is
\f$y_{n+1} = (1-\alpha \Delta t) y_n\f$.  We require for stability that
\f$|y_{n+1}|\le |y_n|\f$, which is to say \f$|1-\alpha \Delta t|\le 1\f$.
Equivalently (since \f$\alpha \Delta t>0\f$),
   \f[\alpha \Delta t\le 2\f]
Therefore we set here
   \f[\Delta t = \frac{2}{\alpha}.\f]
 */
PetscErrorCode PSForceThickness::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {
  PetscErrorCode ierr;
  PetscReal max_dt = 2.0 / alpha * secpera; // convert to seconds

  ierr = input_model->max_timestep(my_t, my_dt, restrict); CHKERRQ(ierr);

  if (restrict) {
    if (max_dt > 0)
      my_dt = PetscMin(max_dt, my_dt);
  }
  else my_dt = max_dt;

  if (my_dt > 0)
    restrict = true;
  else
    restrict = false;

  return 0;
}

//! Adds variables to output files.
void PSForceThickness::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (input_model != NULL)
    input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result["ice_surface_temp"] = ice_surface_temp;
    result["climatic_mass_balance"] = climatic_mass_balance;
  }

  result["ftt_mask"] = ftt_mask.get_metadata();
  result["ftt_target_thk"] = target_thickness.get_metadata();
}

PetscErrorCode PSForceThickness::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ftt_mask")) {
    ierr = ftt_mask.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "ftt_target_thk")) {
    ierr = target_thickness.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype, true); CHKERRQ(ierr);
  }

  ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSForceThickness::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ftt_mask")) {
    ierr = ftt_mask.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "ftt_target_thk")) {
    ierr = target_thickness.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(ice_surface_temp, 0); CHKERRQ(ierr);

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "climatic_mass_balance", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(climatic_mass_balance, 0); CHKERRQ(ierr);

    ierr = ice_surface_mass_flux(tmp); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);

    vars.erase("climatic_mass_balance");
  }

  ierr = input_model->write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}
