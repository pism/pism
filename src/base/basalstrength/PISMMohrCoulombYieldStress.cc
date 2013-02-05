// Copyright (C) 2004--2013 PISM Authors
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

#include "PISMMohrCoulombYieldStress.hh"
#include "PISMHydrology.hh"
#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"

//! \file PISMYieldStress.cc  Process model which computes pseudo-plastic yield stress for the subglacial layer.
/*! \file PISMYieldStress.cc
The output variable of this submodel is \c tauc, the pseudo-plastic yield stress
field that is used in the ShallowStressBalance objects.

In the default implementation PISMMohrCoulombYieldStress [\ref BBssasliding], the
"dry" strength of the (notional) till is a state variable, but private to the
submodel, namely \c tillphi.

Its initialization is nontrivial: either \c -topg_to_phi  heuristic or inverse
modeling so that \c tillphi can be read-in at the beginning of the run. Currently
\c tillphi does not evolve during the run.

This submodel uses a pointer to a PISMHydrology instance to get the basal water
pressure.  Then the effective pressure is combined with tillphi
to compute an updated \c tauc by the Mohr-Coulomb criterion.

This submodel is inactive in floating areas.
*/

PetscErrorCode PISMMohrCoulombYieldStress::allocate() {
  PetscErrorCode ierr;

  ierr = till_phi.create(grid, "tillphi", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = till_phi.set_attrs("model_state",
                            "friction angle for till under grounded ice sheet",
                            "degrees", ""); CHKERRQ(ierr);
  till_phi.time_independent = true; // in this model; need not be
                                    // time-independent in general

  ierr = tauc.create(grid, "tauc", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = tauc.set_attrs("diagnostic",
                        "yield stress for basal till (plastic or pseudo-plastic model)",
                        "Pa", ""); CHKERRQ(ierr);

  ierr = bwp.create(grid, "bwp-in-PISMYieldStress",
                    true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = bwp.set_attrs("internal",
                       "copy of basal water pressure held by PISMMohrCoulombYieldStress",
                       "Pa", ""); CHKERRQ(ierr);
  ierr = Po.create(grid, "overburden_pressure-in-PISMYieldStress",
                    true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = Po.set_attrs("internal",
                       "copy of overburden pressure held by PISMMohrCoulombYieldStress",
                       "Pa", ""); CHKERRQ(ierr);
  return 0;
}


//! Initialize the pseudo-plastic till mechanical model.
/*!
The pseudo-plastic till basal resistance model is governed by this power law
equation,
    \f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}, \f]
where \f$\tau_b=(\tau_{(b)x},\tau_{(b)y})\f$ is the basal shear stress and
\f$U=(u,v)\f$ is the sliding velocity.

We call the scalar field \f$\tau_c(t,x,y)\f$ the \e yield \e stress even when
the power \f$q\f$ is not zero; when that power is zero the formula describes
a plastic material with an actual yield stress.  The constant
\f$U_{\mathtt{th}}\f$ is the \e threshhold \e speed, and \f$q\f$ is the \e pseudo
\e plasticity \e exponent.  The current class computes this yield stress field.
See also IceBasalResistancePlasticLaw::drag().

The strength of the saturated till material, the yield stress, is modeled by a
Mohr-Coulomb relation [\ref Paterson, \ref SchoofStream],
    \f[   \tau_c = c_0 + (\tan \varphi) N, \f]
where \f$N = P_o - P\f$ is the effective pressure on the till (see
PISMHydrology::basal_water_pressure()),

The determination of the till friction angle \f$\varphi(x,y)\f$  is important.
It is assumed in this default model to be a
time-independent factor which describes the strength of the unsaturated "dry"
till material.  Thus it is assumed to change more slowly than the basal water
pressure, and it follows that it changes more slowly than the yield stress and
the basal shear stress.

Option \c -topg_to_phi causes call to topg_to_phi() at the beginning of the run.
This determines the map of \f$\varphi(x,y)\f$.  If this option is note given,
the current method leaves \c tillphi unchanged, and thus either in its
read-in-from-file state or with a default constant value from the config file.
*/
PetscErrorCode PISMMohrCoulombYieldStress::init(PISMVars &vars)
{
  PetscErrorCode ierr;
  bool topg_to_phi_set, plastic_phi_set, bootstrap, i_set,
    tauc_to_phi_set;
  string filename;
  int start;

  variables = &vars;

  ierr = verbPrintf(2, grid.com, "* Initializing the default basal yield stress model...\n"); CHKERRQ(ierr);

  bed_topography = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed_topography == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  ierr = PetscOptionsBegin(grid.com, "", "Options controlling the basal till yield stress model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-plastic_phi", plastic_phi_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-topg_to_phi",
                            "Use the till friction angle parameterization", topg_to_phi_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-tauc_to_phi", "Compute tillphi as a function of tauc and the rest of the model state",
                            tauc_to_phi_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // Get the till friction angle from the the context and ignore options that
  // would be used to set it otherwise.
  IceModelVec2S *till_phi_input = dynamic_cast<IceModelVec2S*>(vars.get("tillphi"));
  if (till_phi_input != NULL) {
    ierr = till_phi.copy_from(*till_phi_input); CHKERRQ(ierr);

    ierr = ignore_option(grid.com, "-plastic_phi"); CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-topg_to_phi"); CHKERRQ(ierr);

    // We do not allow re-gridding in this case.

    return 0;
  }

  if (topg_to_phi_set && plastic_phi_set) {
    PetscPrintf(grid.com, "ERROR: only one of -plastic_phi and -topg_to_phi is allowed.\n");
    PISMEnd();
  }

  if (plastic_phi_set) {

    ierr = till_phi.set(config.get("default_till_phi")); CHKERRQ(ierr);

  } else if (topg_to_phi_set) {

    ierr = verbPrintf(2, grid.com,
                      "  option -topg_to_phi seen; creating tillphi map from bed elev ...\n");
    CHKERRQ(ierr);

    if (i_set || bootstrap) {
      ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);

      PIO nc(grid.com, grid.rank, "guess_mode");
      bool tillphi_present;

      ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
      ierr = nc.inq_var(till_phi.string_attr("short_name"), tillphi_present); CHKERRQ(ierr);
      ierr = nc.close(); CHKERRQ(ierr);

      if (tillphi_present) {
        ierr = verbPrintf(2, grid.com,
                          "PISM WARNING: -topg_to_phi computation will override the '%s' field\n"
                          "              present in the input file '%s'!\n",
                          till_phi.string_attr("short_name").c_str(), filename.c_str()); CHKERRQ(ierr);
      }
    }

    // note option -topg_to_phi will be read again to get comma separated array of parameters
    ierr = topg_to_phi(); CHKERRQ(ierr);

  } else if (i_set || bootstrap) {

    ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);

    if (i_set) {
      ierr = till_phi.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = till_phi.regrid(filename,
                             config.get("bootstrapping_tillphi_value_no_var")); CHKERRQ(ierr);
    }
  }

  ierr = regrid(); CHKERRQ(ierr);

  if (tauc_to_phi_set) {
    string tauc_to_phi_file;
    bool flag;
    ierr = PISMOptionsString("-tauc_to_phi", "Specifies the file tauc will be read from",
                             tauc_to_phi_file, flag, true); CHKERRQ(ierr);

    if (tauc_to_phi_file.empty() == false) {
      // "-tauc_to_phi filename.nc" is given
      ierr = tauc.regrid(tauc_to_phi_file, true); CHKERRQ(ierr);
    } else {
      // "-tauc_to_phi" is given (without a file name); assume that tauc has to
      // be present in an input file
      ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);

      if (bootstrap == false) {
        ierr = tauc.read(filename, start); CHKERRQ(ierr);
      } else {
        ierr = tauc.regrid(filename, true); CHKERRQ(ierr);
      }
    }

    // At this point tauc is initialized in one of the following ways:
    // - from a -tauc_to_phi file
    // - from an input (or bootstrapping) file if -tauc_to_phi did not have an
    //  argument
    //
    // In addition to this, till_phi is initialized
    // - from an input file or
    // - using the -plastic_phi option
    // - using the topg_to_phi option
    //
    // Now tauc_to_phi() will correct till_phi at all locations where grounded
    // ice is present:

    ierr = verbPrintf(2, grid.com, "  Computing till friction angle (tillphi) as a function of the yield stress (tauc)...\n"); 
    CHKERRQ(ierr);

    ierr = tauc_to_phi(); CHKERRQ(ierr);

  } else {
    ierr = tauc.set(0.0); CHKERRQ(ierr);
  }

  // ensure that update() computes tauc at the beginning of the run:
  t = dt = GSL_NAN;

  return 0;
}


PetscErrorCode PISMMohrCoulombYieldStress::regrid() {
  PetscErrorCode ierr;
  bool regrid_file_set, regrid_vars_set;
  string regrid_file;
  set<string> regrid_vars;

  ierr = PetscOptionsBegin(grid.com, "", "PISMMohrCoulombYieldStress regridding options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "regridding file name",
                             regrid_file, regrid_file_set); CHKERRQ(ierr);
    ierr = PISMOptionsStringSet("-regrid_vars", "comma-separated list of regridding variables",
                                "", regrid_vars, regrid_vars_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (! regrid_file_set) return 0;

  // stop if the user did not ask to regrid tillphi
  if (!set_contains(regrid_vars, till_phi.string_attr("short_name")))
    return 0;

  ierr = till_phi.regrid(regrid_file, true); CHKERRQ(ierr);

  return 0;
}



void PISMMohrCoulombYieldStress::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["tillphi"] = till_phi.get_metadata();
}

PetscErrorCode PISMMohrCoulombYieldStress::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  if (set_contains(vars, "tillphi")) {
    PetscErrorCode ierr = till_phi.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMMohrCoulombYieldStress::write_variables(set<string> vars, const PIO &nc) {
  if (set_contains(vars, "tillphi")) {
    PetscErrorCode ierr = till_phi.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


//! Update the till yield stress for use in the pseudo-plastic till basal stress
//! model.  See also IceBasalResistancePlasticLaw.
/*!
Updates yield stress \f$\tau_c\f$ based on modeled basal water pressure.  We implement
formula (2.4) in [\ref SchoofStream], the Mohr-Coulomb criterion:
    \f[   \tau_c = \mu (P_o - P), \f]
where \f$\tau_c\f$ is the till yield stress, \f$P_o\f$ is the ice over-burden
pressure, \f$P\f$ is the modeled basal water pressure, and \f$\mu\f$ is a
strength coefficient for the mineral till (at least, it is independent of
\f$p_w\f$).  The difference
    \f[   N = P_o - P   \f]
is the effective pressure on the till.  Both \f$P_o,P\f$ are provided by
PISMHydrology.

We modify Schoof's formula by allowing a small till cohesion \f$c_0\f$
and by expressing the coefficient as the tangent of a till friction angle
\f$\varphi\f$:
    \f[   \tau_c = c_0 + (\tan \varphi) N. \f]
Option  \c -plastic_c0 controls it \f$c_0\f$; see [\ref Paterson] table 8.1
regarding values.
 */
PetscErrorCode PISMMohrCoulombYieldStress::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t = my_t; dt = my_dt;
  // this model performs a "diagnostic" computation (i.e. without time-stepping)

  bool use_ssa_when_grounded = config.get_flag("use_ssa_when_grounded");

  const PetscScalar
    high_tauc = config.get("high_tauc");

  if (hydrology) {
    ierr = hydrology->subglacial_water_pressure(bwp); CHKERRQ(ierr);
    ierr = hydrology->overburden_pressure(Po); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 3,
            "PISM ERROR: PISMHydrology* hydrology is NULL in PISMMohrCoulombYieldStress::update()");
  }

  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = tauc.begin_access(); CHKERRQ(ierr);
  ierr = bwp.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  ierr = till_phi.begin_access(); CHKERRQ(ierr);

  MaskQuery m(*mask);

  PetscInt GHOSTS = grid.max_stencil_width;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      if (use_ssa_when_grounded == false) {
        if (m.grounded(i, j)) {
          // large yield stress if grounded and -ssa_floating_only is set
          tauc(i, j) = high_tauc;
        } else {
          tauc(i, j) = 0.0;
        }
        continue;
      }
      if (m.ocean(i, j)) {
        tauc(i, j) = 0.0;
      } else if (m.ice_free(i, j)) {
        tauc(i, j) = high_tauc;  // large yield stress if grounded and ice-free
      } else { // grounded and there is some ice
        const PetscScalar N = Po(i,j) - bwp(i,j);
        tauc(i, j) = till_c_0 + N * tan((pi/180.0) * till_phi(i, j));
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc.end_access(); CHKERRQ(ierr);
  ierr = till_phi.end_access(); CHKERRQ(ierr);
  ierr = bwp.end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMMohrCoulombYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  return tauc.copy_to(result);
}


//! Computes the till friction angle phi as a piecewise linear function of bed elevation, according to user options.
/*!
Computes the till friction angle \f$\phi(x,y)\f$ at a location, namely
\c IceModel::vtillphi, as the following increasing, piecewise-linear function of
the bed elevation \f$b(x,y)\f$.  Let
	\f[ M = (\phi_{\text{max}} - \phi_{\text{min}}) / (b_{\text{max}} - b_{\text{min}}) \f]
be the slope of the nontrivial part.  Then
	\f[ \phi(x,y) = \begin{cases}
	        \phi_{\text{min}}, & b(x,y) \le b_{\text{min}}, \\
	        \phi_{\text{min}} + (b(x,y) - b_{\text{min}}) \,M,
	                          &  b_{\text{min}} < b(x,y) < b_{\text{max}}, \\
	        \phi_{\text{max}}, & b_{\text{max}} \le b(x,y), \end{cases} \f]
The exception is if the point is marked as floating, in which case the till friction angle
is set to the value \c phi_ocean.

The default values are vaguely suitable for Antarctica, perhaps:
- \c phi_min = 5.0 degrees,
- \c phi_max = 15.0 degrees,
- \c topg_min = -1000.0 m,
- \c topg_max = 1000.0 m,
- \c phi_ocean = 10.0 degrees.

If the user gives option <code>-topg_to_phi A,B,C,D</code> then \c phi_ocean
is not used. Instead, the same rule as above for grounded ice is used.
 */
PetscErrorCode PISMMohrCoulombYieldStress::topg_to_phi() {
  PetscErrorCode ierr;
  bool  topg_to_phi_set;
  vector<double> inarray(4);

  // default values:
  inarray[0] = 5.0;
  inarray[1] = 15.0;
  inarray[2] = -1000.0;
  inarray[3] = 1000.0;

  // read the comma-separated list of four values
  ierr = PISMOptionsRealArray("-topg_to_phi", "phi_min, phi_max, topg_min, topg_max",
                              inarray, topg_to_phi_set); CHKERRQ(ierr);

  if (topg_to_phi_set == false) {
    SETERRQ(grid.com, 1, "HOW DID I GET HERE? ... ending...\n");
  }

  if (inarray.size() != 4) {
    PetscPrintf(grid.com,
                "PISM ERROR: option -topg_to_phi requires a comma-separated list with 4 numbers; got %d\n",
                inarray.size());
    PISMEnd();
  }

  PetscReal phi_min = inarray[0], phi_max = inarray[1],
    topg_min = inarray[2], topg_max = inarray[3];

  if (phi_min >= phi_max) {
    PetscPrintf(grid.com,
                "PISM ERROR: invalid -topg_to_phi arguments: phi_min < phi_max is required\n");
    PISMEnd();
  }

  if (topg_min >= topg_max) {
    PetscPrintf(grid.com,
                "PISM ERROR: invalid -topg_to_phi arguments: topg_min < topg_max is required\n");
    PISMEnd();
  }

  ierr = verbPrintf(2, grid.com,
                    "  till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
                    "            /  %5.2f                                 for   topg < %.f\n"
                    "      phi = |  %5.2f + (topg - (%.f)) * (%.2f / %.f)   for   %.f < topg < %.f\n"
                    "            \\  %5.2f                                 for   %.f < topg\n",
                    phi_min, topg_min,
                    phi_min, topg_min, phi_max - phi_min, topg_max - topg_min, topg_min, topg_max,
                    phi_max, topg_max); CHKERRQ(ierr);

  PetscReal slope = (phi_max - phi_min) / (topg_max - topg_min);

  ierr = bed_topography->begin_access(); CHKERRQ(ierr);
  ierr = till_phi.begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      PetscScalar bed = (*bed_topography)(i, j);

      if (bed <= topg_min) {
        till_phi(i, j) = phi_min;
      } else if (bed >= topg_max) {
        till_phi(i, j) = phi_max;
      } else {
        till_phi(i, j) = phi_min + (bed - topg_min) * slope;
      }

    }
  }

  ierr = bed_topography->end_access(); CHKERRQ(ierr);
  ierr = till_phi.end_access(); CHKERRQ(ierr);

  // communicate ghosts so that the tauc computation can be performed locally
  // (including ghosts of tauc, that is)
  ierr = till_phi.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMMohrCoulombYieldStress::tauc_to_phi() {
  PetscErrorCode ierr;

  if (hydrology) {
    ierr = hydrology->subglacial_water_pressure(bwp); CHKERRQ(ierr);
    ierr = hydrology->overburden_pressure(Po); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 3,
            "PISM ERROR: PISMHydrology* hydrology is NULL in PISMMohrCoulombYieldStress::tauc_to_phi()");
  }

  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = tauc.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  ierr = bwp.begin_access(); CHKERRQ(ierr);
  ierr = till_phi.begin_access(); CHKERRQ(ierr);

  MaskQuery m(*mask);

  PetscInt GHOSTS = grid.max_stencil_width;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {

      if (m.ocean(i, j)) {
        // no change
      } else if (m.ice_free(i, j)) {
        // no change
      } else { // grounded and there is some ice
        const PetscScalar N = Po(i,j) - bwp(i,j);
        till_phi(i, j) = 180.0/pi * atan((tauc(i, j) - till_c_0) / N);
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc.end_access(); CHKERRQ(ierr);
  ierr = till_phi.end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  ierr = bwp.end_access(); CHKERRQ(ierr);

  return 0;
}
