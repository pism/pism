// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "PSForceThickness.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "PIO.hh"
#include "PISMConfig.hh"
#include "Mask.hh"

#include "pism_options.hh"

#include "error_handling.hh"

namespace pism {

///// "Force-to-thickness" mechanism
PSForceThickness::PSForceThickness(IceGrid &g, const Config &conf, SurfaceModel *input)
  : PSModifier(g, conf, input),
    m_climatic_mass_balance(g.get_unit_system(), "climatic_mass_balance", grid),
    m_climatic_mass_balance_original(g.get_unit_system(), "climatic_mass_balance_original", grid),
    m_ice_surface_temp(g.get_unit_system(), "ice_surface_temp", grid) {

  m_ice_thickness = NULL;
  m_alpha = config.get("force_to_thickness_alpha", "yr-1", "s-1");
  m_alpha_ice_free_factor = config.get("force_to_thickness_ice_free_alpha_factor");
  m_ice_free_thickness_threshold = config.get("force_to_thickness_ice_free_thickness_threshold");

  m_target_thickness.create(grid, "thk", WITHOUT_GHOSTS);
  // will set attributes in init()

  m_ftt_mask.create(grid, "ftt_mask", WITHOUT_GHOSTS);
  m_ftt_mask.set_attrs("diagnostic",
                       "mask specifying where to apply the force-to-thickness mechanism",
                       "", ""); // no units and no standard name
  m_ftt_mask.set(1.0); // default: applied in whole domain
  m_ftt_mask.write_in_glaciological_units = true;

  m_climatic_mass_balance.set_string("pism_intent", "diagnostic");
  m_climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  m_climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.set_units("kg m-2 s-1");
  m_climatic_mass_balance.set_glaciological_units("kg m-2 year-1");

  m_climatic_mass_balance_original.set_string("pism_intent", "diagnostic");
  m_climatic_mass_balance_original.set_string("long_name",
                                            "surface mass balance rate before the adjustment using -surface ...,forcing");
  m_climatic_mass_balance_original.set_units("kg m-2 s-1");
  m_climatic_mass_balance_original.set_glaciological_units("kg m-2 year-1");

  m_ice_surface_temp.set_string("pism_intent", "diagnostic");
  m_ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  m_ice_surface_temp.set_units("K");
}

PSForceThickness::~PSForceThickness() {
  // empty
}

void PSForceThickness::attach_atmosphere_model(AtmosphereModel *input) {
  input_model->attach_atmosphere_model(input);
}

void PSForceThickness::init(Vars &vars) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  input_model->init(vars);

  verbPrintf(2, grid.com,
             "* Initializing force-to-thickness mass-balance modifier...\n");

  bool file_set = false;
  OptionsString("-force_to_thickness_file",
                "Specifies the target thickness file for the force-to-thickness mechanism",
                m_input_file, file_set, false);

  if (file_set == false) {
    throw RuntimeError("surface model forcing requires the -force_to_thickness_file option.");
  }

  double ftt_alpha = grid.convert(m_alpha, "s-1", "yr-1");
  bool ftt_alpha_set = false;
  OptionsReal("-force_to_thickness_alpha",
              "Specifies the value of force-to-thickness alpha in per-year units",
              ftt_alpha, ftt_alpha_set);

  bool ftt_alpha_ice_free_set = false;
  OptionsReal("-force_to_thickness_ice_free_alpha_factor",
              "Set the multiplicative factor for alpha to use in ice-free areas",
              m_alpha_ice_free_factor, ftt_alpha_ice_free_set);

  bool ftt_ice_free_thickness_threshold_set = false;
  OptionsReal("-force_to_thickness_ice_free_thickness_threshold",
              "Specifies the ice thickness threshold used to determine whether a location is ice-free, in m",
              m_ice_free_thickness_threshold, ftt_ice_free_thickness_threshold_set);

  m_ice_thickness = vars.get_2d_scalar("land_ice_thickness");
  m_pism_mask     = vars.get_2d_mask("mask");

  // determine exponential rate alpha from user option or from factor; option
  // is given in a^{-1}
  if (ftt_alpha_set == true) {
    verbPrintf(3, grid.com, "    option -force_to_thickness_alpha seen\n");
    m_alpha = grid.convert(ftt_alpha, "yr-1", "s-1");
  }

  verbPrintf(2, grid.com,
             "    alpha = %.6f year-1 for -force_to_thickness mechanism\n"
             "    alpha = %.6f year-1 in areas with target ice thickness of less than %.3f meters\n",
             grid.convert(m_alpha, "s-1", "yr-1"),
             m_alpha_ice_free_factor * grid.convert(m_alpha, "s-1", "yr-1"),
             m_ice_free_thickness_threshold);

  // m_input_file now contains name of -force_to_thickness file; now check
  // it is really there; and regrid the target thickness
  PIO nc(grid, "guess_mode");
  bool mask_exists = false;
  nc.open(m_input_file, PISM_READONLY);
  mask_exists = nc.inq_var("ftt_mask");
  nc.close();

  verbPrintf(2, grid.com,
             "    reading target thickness 'thk' from %s ...\n"
             "    (this field will appear in output file as 'ftt_target_thk')\n",
             m_input_file.c_str());
  {
    m_target_thickness.set_name("thk"); // name to read by
    // set attributes for the read stage; see below for reset
    m_target_thickness.set_attrs("diagnostic",
                                 "target thickness for force-to-thickness mechanism (hit this at end of run)",
                                 "m",
                                 "land_ice_thickness"); // standard_name *to read by*

    m_target_thickness.regrid(m_input_file, CRITICAL);

    // reset name to avoid confusion; set attributes again to overwrite "read by" choices above
    m_target_thickness.set_name("ftt_target_thk");
    m_target_thickness.set_attrs("diagnostic",
                                 "target thickness for force-to-thickness mechanism (wants to hit this at end of run)",
                                 "m",
                                 "");  // no CF standard_name, to put it mildly

    m_target_thickness.write_in_glaciological_units = true;
  }

  if (mask_exists) {
    verbPrintf(2, grid.com,
               "    reading force-to-thickness mask 'ftt_mask' from %s ...\n",
               m_input_file.c_str());
    m_ftt_mask.regrid(m_input_file, CRITICAL);
  }
}

/*!
If `-force_to_thickness_file` `foo.nc` is in use then vthktarget will have a target ice thickness
map.  Let \f$H_{\text{tar}}\f$ be this target thickness,
and let \f$H\f$ be the current model thickness.  Recall that the mass continuity
equation solved by IceModel::massContExplicitStep() is
  \f[ \frac{\partial H}{\partial t} = M - S - \nabla\cdot \mathbf{q} \f]
and that this procedure is supposed to produce \f$M\f$.
In this context, the semantics of `-force_to_thickness` are that \f$M\f$ is modified
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

The constant \f$\alpha\f$ has a default value `pism_config:force_to_thickness_alpha`.

The next example uses files generated from the EISMINT-Greenland experiment;
see the corresponding chapter of the User's Manual.

Suppose we regard the SSL2 run as a spin-up to reach a better temperature field.
It is a spinup in which the surface was allowed to evolve.  Assume the
early file `green20km_y1.nc` has the target thickness, because it essentially
has the input thickness.  This script adds a 500 a run, to finalize the spinup,
in which the ice sheet geometry goes from the the thickness values in
`green_ssl2_110ka.nc` to values very close to those in `green20km_y1.nc`:
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
#  $ export PISM_EXEC="pismr -energy cold"
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
  -surface pdd,forcing -pdd_fausto -force_to_thickness_file green20km_y1.nc \
  -o default_force.nc -ts_file ts_default_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -1000.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd,forcing -pdd_fausto -force_to_thickness_file green20km_y1.nc -force_to_thickness_alpha 0.005 \
    -o weak_force.nc -ts_file ts_weak_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd


cmd="$PISM_MPIDO $NN $PISM -ys -1000.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd,forcing -pdd_fausto -force_to_thickness_file green20km_y1.nc -force_to_thickness_alpha 0.05 \
    -o strong_force.nc -ts_file ts_strong_force.nc -ts_times -1000:yearly:0"
$PISM_DO $cmd

\endcode
The script also has a run with no forcing, one with forcing at a lower alpha value,
a factor of five smaller than the default, and one with a forcing at a higher alpha value, a factor of five higher.
 */
void PSForceThickness::ice_surface_mass_flux(IceModelVec2S &result) {

  // get the surface mass balance result from the next level up
  input_model->ice_surface_mass_flux(result);

  verbPrintf(5, grid.com,
             "    updating surface mass balance using -force_to_thickness mechanism ...");

  double ice_density = config.get("ice_density");

  MaskQuery m(*m_pism_mask);

  IceModelVec::AccessList list;
  list.add(*m_pism_mask);
  list.add(*m_ice_thickness);
  list.add(m_target_thickness);
  list.add(m_ftt_mask);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_ftt_mask(i,j) > 0.5 && m.grounded(i, j)) {
      if (m_target_thickness(i,j) >= m_ice_free_thickness_threshold) {
        result(i,j) += ice_density * m_alpha * (m_target_thickness(i,j) - (*m_ice_thickness)(i,j));
      } else {
        result(i,j) += ice_density * m_alpha * m_alpha_ice_free_factor * (m_target_thickness(i,j) - (*m_ice_thickness)(i,j));
      }
    }
  }
  // no communication needed
}

//! Does not modify ice surface temperature.
void PSForceThickness::ice_surface_temperature(IceModelVec2S &result) {
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
void PSForceThickness::max_timestep(double my_t, double &my_dt, bool &restrict) {
  double max_dt = grid.convert(2.0 / m_alpha, "years", "seconds");

  input_model->max_timestep(my_t, my_dt, restrict);

  if (restrict) {
    if (max_dt > 0) {
      my_dt = std::min(max_dt, my_dt);
    }
  } else {
    my_dt = max_dt;
  }

  if (my_dt > 0) {
    restrict = true;
  } else {
    restrict = false;
  }
}

//! Adds variables to output files.
void PSForceThickness::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (input_model != NULL) {
    input_model->add_vars_to_output(keyword, result);
  }

  if (keyword == "medium" || keyword == "big") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
    result.insert("climatic_mass_balance_original");
  }

  result.insert("ftt_mask");
  result.insert("ftt_target_thk");
}

void PSForceThickness::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  if (set_contains(vars, "ftt_mask")) {
    m_ftt_mask.define(nc, nctype);
  }

  if (set_contains(vars, "ftt_target_thk")) {
    m_target_thickness.define(nc, nctype);
  }

  if (set_contains(vars, "ice_surface_temp")) {
    m_ice_surface_temp.define(nc, nctype, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    m_climatic_mass_balance.define(nc, nctype, true);
  }

  if (set_contains(vars, "climatic_mass_balance_original")) {
    m_climatic_mass_balance_original.define(nc, nctype, true);
  }

  input_model->define_variables(vars, nc, nctype);
}

void PSForceThickness::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ftt_mask")) {
    m_ftt_mask.write(nc);
  }

  if (set_contains(vars, "ftt_target_thk")) {
    m_target_thickness.write(nc);
  }

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = m_ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);

    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance_original")) {
    IceModelVec2S tmp;
    tmp.create(grid, "climatic_mass_balance_original", WITHOUT_GHOSTS);
    tmp.metadata() = m_climatic_mass_balance_original;

    input_model->ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("climatic_mass_balance_original");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    tmp.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS);
    tmp.metadata() = m_climatic_mass_balance;

    ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("climatic_mass_balance");
  }

  input_model->write_variables(vars, nc);
}

} // end of namespace pism
