// Copyright (C) 2004--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "PISMYieldStress.hh"
#include "Mask.hh"

/*! \file PISMYieldStress.cc
The output variable of this submodel is tauc.  It is used in the SSA objects.

The "dry" strength of the (notional) till is a state variable, but private to
the submodel: tillphi. Its initialization is nontrivial: either -topg_to_phi
heuristic or inverse modeling (to be implemented ...). Currently tillphi does
not evolve during the run.

This submodel uses vHmelt as an input at each update. A basal water pressure is
computed from vHmelt, then that basal water pressure is combined with tillphi
to compute an updated tauc.

This submodel is inactive in floating areas.
*/

PetscErrorCode PISMDefaultYieldStress::allocate() {
  PetscErrorCode ierr;

  ierr = till_phi.create(grid, "tillphi", true); CHKERRQ(ierr);
  ierr = till_phi.set_attrs("model_state",
                            "friction angle for till under grounded ice sheet",
                            "degrees", ""); CHKERRQ(ierr);
  till_phi.time_independent = true; // in this model; need not be
                                    // time-independent in general

  return 0;
}

//! Initialize the pseudo-plastic till mechanical model.
/*! 
The pseudo-plastic till basal resistance model is governed by this power law
equation,
    \f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}, \f]
where \f$\tau_b=(\tau_{(b)x},\tau_{(b)y})\f$ is the basal shear stress and
\f$U=(u,v)\f$ is the sliding velocity.  We call the scalar field 
\f$\tau_c(t,x,y)\f$ the \e pseudo \e yield \e stress.  The constant
\f$U_{\mathtt{th}}\f$ is the \e threshhold \e speed, and \f$q\f$ is the \e pseudo
\e plasticity \e exponent.  See IceBasalResistancePlasticLaw::drag().  See also
updateYieldStressUsingBasalWater() and getBasalWaterPressure() and diffuseHmelt()
for important model equations.  

Because the strength of the saturated till material is modeled by a Mohr-Coulomb
relation [\ref Paterson, \ref SchoofStream],
    \f[   \tau_c = c_0 + (\tan \varphi) N, \f]
where \f$N = \rho g H - p_w\f$ is the effective pressure on the till (see
getBasalWaterPressure()), the determination of the till friction angle
\f$\varphi(x,y)\f$  is important.  The till friction angle is assumed to be a
time-independent factor which describes the strength of the unsaturated "dry"
till material.  Thus it is assumed to change more slowly than the basal water 
pressure.  It follows that it changes more slowly than the yield stress and/or
the basal shear stress.

The current method determines the map of till friction angle \c vtillphi
according to options.  Option \c -topg_to_phi causes call to
topg_to_phi().  If neither option is given, the current method
leaves \c tillphi unchanged, and thus either in its read-in-from-file state or 
with a default constant value from the config file.

The current method reads the experimental parameter option \c -sliding_scale.
A scale factor of \f$A\f$ is intended to increase basal sliding rate by
\f$A\f$.  It would have exactly this effect \e if the driving stress were
\e hypothetically completely held by the basal resistance.  Thus this scale factor
is used to reduce (if \c -sliding_scale \f$A\f$ with \f$A > 1\f$) or increase
(if \f$A < 1\f$) the value of the (pseudo-) yield stress \c tauc.  The concept
behind this is described at
http://websrv.cs.umt.edu/isis/index.php/Category_1:_Whole_Ice_Sheet#Initial_Experiment_-_E1_-_Increased_Basal_Lubrication.
Specifically, the concept behind this mechanism is to suppose equality of driving
and basal shear stresses,
    \f[ \rho g H \nabla h = \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}. \f]
(<i>For emphasis:</i> The membrane stress held by the ice itself is missing from
this incomplete stress balance.)  Thus the pseudo yield stress
\f$\tau_c\f$ would be related to the sliding speed \f$|\mathbf{U}|\f$ by
  \f[ |\mathbf{U}| = \frac{C}{\tau_c^{1/q}} \f]
for some (geometry-dependent) constant \f$C\f$.  Multiplying \f$|\mathbf{U}|\f$
by \f$A\f$ in this equation corresponds to dividing \f$\tau_c\f$ by \f$A^q\f$.
The current method sets-up the mechanism, and updateYieldStressUsingBasalWater()
actually computes it.  Note that the mechanism has no effect whatsoever if 
\f$q=0\f$, which is the purely plastic case. In that case there is \e no direct
relation between the yield stress and the sliding velocity, and the difference 
between the driving stress and the yield stress is entirely held by the membrane
stresses.  (There is also no singular mathematical operation as \f$A^q = A^0 = 1\f$.)
*/

PetscErrorCode PISMDefaultYieldStress::init(PISMVars &vars) {
  PetscErrorCode ierr;
  PetscScalar pseudo_plastic_q = config.get("pseudo_plastic_q");
  bool topg_to_phi_set, plastic_phi_set, bootstrap, i_set;
  string filename;
  int start;

  variables = &vars;

  ierr = till_phi.set(config.get("default_till_phi")); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Options controlling the basal till yield stress model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-plastic_phi", plastic_phi_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-topg_to_phi",
                            "Use the till friction angle parameterization", topg_to_phi_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);

    bool scaleSet;
    double slidescale;
    ierr = PISMOptionsReal("-sliding_scale",
                           "Divides pseudo-plastic tauc (yield stress) by given factor;"
                           " this would increase sliding by given factor in absence of membrane stresses",
                           slidescale, scaleSet); CHKERRQ(ierr);
    if (scaleSet) { // only modify config if option set; otherwise leave alone
      if (slidescale > 0.0) {
        ierr = verbPrintf(2, grid.com,
                          "option -sliding_scale read; pseudo yield stress tauc will be divided by %.4f to\n"
                          "  cause notional sliding speed-up by given factor %.4f ...\n",
                          pow(slidescale, pseudo_plastic_q), slidescale);  CHKERRQ(ierr);
        sliding_scale = slidescale;
      } else {
        ierr = verbPrintf(1, grid.com,
                          "PISM WARNING: negative or zero value given for option -sliding_scale ignored\n");
        CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (i_set || bootstrap) {
    ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);

    if (i_set) {
      ierr = till_phi.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = till_phi.regrid(filename,
                             config.get("bootstrapping_tillphi_value_no_var")); CHKERRQ(ierr);
    }
  }

  if (topg_to_phi_set && plastic_phi_set) {
    PetscPrintf(grid.com, "ERROR: only one of -plastic_phi and -topg_to_phi is allowed.\n");
    PISMEnd();
  }

  basal_water_thickness = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (basal_water_thickness == NULL) SETERRQ(1, "bwat is not available");

  basal_melt_rate = dynamic_cast<IceModelVec2S*>(vars.get("bmelt"));
  if (basal_melt_rate == NULL) SETERRQ(1, "bmelt is not available");

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (ice_thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  bed_topography = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed_topography == NULL) SETERRQ(1, "bedrock_altitude is not available");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");

  if (topg_to_phi_set) {
    ierr = verbPrintf(2, grid.com,
                      "option -topg_to_phi seen; creating till friction angle map from bed elevation...\n");
    CHKERRQ(ierr);
    // note option -topg_to_phi will be read again to get comma separated array of parameters
    ierr = topg_to_phi(); CHKERRQ(ierr);
  }

  ierr = regrid(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMDefaultYieldStress::regrid() {
  PetscErrorCode ierr;
  bool regrid_file_set, regrid_vars_set;
  string regrid_file;
  vector<string> regrid_vars;

  ierr = PetscOptionsBegin(grid.com, "", "PISMDefaultYieldStress regridding options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "regridding file name",
                             regrid_file, regrid_file_set); CHKERRQ(ierr);
    ierr = PISMOptionsStringArray("-regrid_vars", "comma-separated list of regridding variables",
                                  "", regrid_vars, regrid_vars_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (! regrid_file_set) return 0;

  set<string> vars;
  for (unsigned int i = 0; i < regrid_vars.size(); ++i)
    vars.insert(regrid_vars[i]);

  // stop if the user did not ask to regrid tillphi
  if (!set_contains(vars, till_phi.string_attr("short_name")))
    return 0;

  ierr = till_phi.regrid(regrid_file, true); CHKERRQ(ierr);

  return 0;
}



void PISMDefaultYieldStress::add_vars_to_output(string /*keyword*/, set<string> &result) {
  result.insert("tillphi");
}


void PISMDefaultYieldStress::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  dict["bwp"] = new PYS_bwp(this, grid, *variables);
}


PetscErrorCode PISMDefaultYieldStress::define_variables(set<string> vars, const NCTool &nc,
                                                 nc_type nctype) {
  if (set_contains(vars, "tillphi")) {
    PetscErrorCode ierr = till_phi.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMDefaultYieldStress::write_variables(set<string> vars, string filename) {
  if (set_contains(vars, "tillphi")) {
    PetscErrorCode ierr = till_phi.write(filename); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode PISMDefaultYieldStress::update(PetscReal t_years, PetscReal dt_years) {
  t = t_years; dt = dt_years;
  // this model performs a "diagnostic" computation (i.e. without time-stepping)
  return 0;
}

//! Update the till yield stress for the pseudo-plastic till SSA model.
/*!
Updates based on modeled basal water pressure.  We implement
formula (2.4) in [\ref SchoofStream],
    \f[   \tau_c = \mu (\rho g H - p_w), \f]
where \f$\tau_c\f$ is the till yield stress, \f$\rho g H\f$ is the ice over-burden
pressure (in the shallow approximation), \f$p_w\f$ is the modeled
pore (=basal) water pressure, and \f$\mu\f$ is a strength coefficient for the
mineral till (at least, it is independent of \f$p_w\f$).  The difference
    \f[   N = \rho g H - p_w   \f]
is the effective pressure on the till.
 
We modify Schoof's formula by allowing a small till cohesion \f$c_0\f$
and by expressing the coefficient as the tangent of a till friction angle
\f$\varphi\f$:
    \f[   \tau_c = c_0 + (\tan \varphi) N. \f]
See [\ref Paterson] table 8.1) regarding values of \f$c_0\f$.
Option  \c -plastic_c0 controls it.

The major concern with this is the model for basal water pressure \f$p_w\f$.  
See getBasalWaterPressure().  See also [\ref BBssasliding] for a discussion
of a complete model using these tools.

Note that IceModel::updateSurfaceElevationAndMask() also
checks whether use_ssa_when_grounded is true and if so it sets all mask points to
DRAGGING.
 */
PetscErrorCode PISMDefaultYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  PetscErrorCode ierr;

  bool use_ssa_when_grounded = config.get_flag("use_ssa_when_grounded");
  // only makes sense when use_ssa_when_grounded == TRUE
  if (use_ssa_when_grounded == PETSC_FALSE) {
    SETERRQ(1, "use_ssa_when_grounded == PETSC_FALSE but updateYieldStressFromHmelt() called");
  }


  PetscScalar till_pw_fraction = config.get("till_pw_fraction"),
    till_c_0 = config.get("till_c_0") * 1e3, // convert from kPa to Pa
    hmelt_max = config.get("hmelt_max"),
    ice_density = config.get("ice_density"),
    standard_gravity = config.get("standard_gravity");


  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = basal_water_thickness->begin_access(); CHKERRQ(ierr);
  ierr = basal_melt_rate->begin_access(); CHKERRQ(ierr);
  ierr = till_phi.begin_access(); CHKERRQ(ierr);

  MaskQuery m(*mask);

  PetscInt GHOSTS = grid.max_stencil_width;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      if (m.ocean(i, j)) {
        result(i, j) = 0.0;
      } else if (m.ice_free(i, j)) {
        result(i, j) = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar if no ice
      } else { // grounded and there is some ice
        const PetscScalar
          p_over = ice_density * standard_gravity * (*ice_thickness)(i, j), // FIXME task #7297
          p_w    = basal_water_pressure((*ice_thickness)(i, j),
                                        (*basal_water_thickness)(i, j),
                                        (*basal_melt_rate)(i, j),
                                        till_pw_fraction, hmelt_max),
          N      = p_over - p_w;  // effective pressure on till

        result(i, j) = till_c_0 + N * tan((pi/180.0) * till_phi(i, j));
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = till_phi.end_access(); CHKERRQ(ierr);
  ierr = basal_melt_rate->end_access(); CHKERRQ(ierr);
  ierr = basal_water_thickness->end_access(); CHKERRQ(ierr);

  // scale tauc if desired
  if (sliding_scale > 0.0) {
    const PetscScalar q = config.get("pseudo_plastic_q");
    result.scale(1.0 / pow(sliding_scale, q));
  }

  return 0;
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
PetscErrorCode PISMDefaultYieldStress::topg_to_phi() {
  PetscErrorCode ierr;

  PetscInt    Nparam=5;
  PetscReal   inarray[5] = {5.0, 15.0, -1000.0, 1000.0, 10.0};

  // read comma-separated array of zero to five values
  PetscTruth  topg_to_phi_set;
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-topg_to_phi", inarray, &Nparam, &topg_to_phi_set);
  CHKERRQ(ierr);
  if (topg_to_phi_set != PETSC_TRUE) {
    SETERRQ(1, "HOW DID I GET HERE? ... ending...\n");
  }

  if ((Nparam > 5) || (Nparam < 4)) {
    ierr = verbPrintf(1, grid.com,
                      "PISM ERROR: option -topg_to_phi provided with more than 5 or fewer than 4\n"
                      "            arguments ... ENDING ...\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  PetscReal phi_min = inarray[0], phi_max = inarray[1],
    topg_min = inarray[2], topg_max = inarray[3],
    phi_ocean = inarray[4];

  ierr = verbPrintf(2, grid.com,
                    "  till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
                    "            /  %5.2f                                 for   topg < %.f\n"
                    "      phi = |  %5.2f + (topg - %.f) * (%.2f / %.f)   for   %.f < topg < %.f\n"
                    "            \\  %5.2f                                 for   %.f < topg\n",
                    phi_min, topg_min,
                    phi_min, topg_min, phi_max - phi_min, topg_max - topg_min, topg_min, topg_max,
                    phi_max, topg_max); CHKERRQ(ierr);

  if (Nparam == 5) {
    ierr = verbPrintf(2, grid.com,
                      "      (also using phi = %5.2f in floating ice or ice free ocean)\n",
                      phi_ocean); CHKERRQ(ierr);
  }

  MaskQuery m(*mask);

  PetscReal slope = (phi_max - phi_min) / (topg_max - topg_min);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = bed_topography->begin_access(); CHKERRQ(ierr);
  ierr = till_phi.begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      PetscScalar bed = (*bed_topography)(i, j);

      if (m.grounded(i, j) || (Nparam < 5)) {
        if (bed <= topg_min) {
          till_phi(i, j) = phi_min;
        } else if (bed >= topg_max) {
          till_phi(i, j) = phi_max;
        } else {
          till_phi(i, j) = phi_min + (bed - topg_min) * slope;
        }
      } else {
        till_phi(i,j) = phi_ocean;
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = bed_topography->end_access(); CHKERRQ(ierr);
  ierr = till_phi.end_access(); CHKERRQ(ierr);

  // communicate ghosts so that the tauc computation can be performed locally
  // (including ghosts of tauc, that is)
  ierr = till_phi.beginGhostComm(); CHKERRQ(ierr);
  ierr = till_phi.endGhostComm(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute modeled pressure in subglacial liquid water using thickness
//! \brief   of water layer, and possibly the melt rate, at the base.
/*!
This procedure provides a simple model of basal water pressure \f$p_w\f$, which
is modeled as a function of the thickness of the basal stored water plus 
(optionally) the basal melt rate.

Input \c bwat is thickness of basal water.  Input \c bmr is the basal melt rate.
Because both \c bwat and \c bmr are zero at points where base of ice is 
below the pressure-melting temperature, the modeled basal water pressure is
zero when the base is frozen.

The inequality \c bwat \f$\le\f$ \c hmelt_max is required at input, and an
error is thrown if not.

Regarding the physics, compare the water pressure computed by formula (4) in
[\ref Ritzetal2001], where the pressure is a function of sea level and bed
elevation.  Also, the method using "elevation of the bed at the grounding line"
as in [\ref LingleBrown1987] is not implementable because that elevation is
at an unknowable location.  (We are not doing a flow line model!)

Several options control the water pressure model:
  - \c -[no_]\c bmr_enhance  toggle the basal melt rate dependency in water
                           pressure; DEFAULT IS OFF
  - \c -bmr_enhance_scale  sets the value; defaults to 0.10, since 10 cm/a is
                           a significant enough level of basal melt rate to
                           cause weakening effect for saturated till;
                           argument is in m/a; must set \c -bmr_enhance for this
                           to have effect
  - \c -plastic_pwfrac     controls parameter till_pw_fraction
  - \c -[no_]\c thk_eff    toggle the thickness effect: for smaller thicknesses there
                           is a reduction in basal water pressure, a conceptual
                           near-margin drainage mechanism; DEFAULT IS OFF
  - \c -thk_eff_reduced    factor by which basal water pressure is reduced; default is
                           0.97; must set \c -thk_eff for this to have effect
  - \c -thk_eff_H_high     maximum thickness at which effect is applied; default
                           is 2000 m; must set \c -thk_eff for this to have
                           any effect
  - \c -thk_eff_H_low      thickness at which thickness effect is full strength;
                           default is 1000 m; must set \c -thk_eff for this to
                           have any effect

At several places in the code the effective pressure on the mineral part of the
till is computed by these lines, which are recommended for this purpose:

<code>
  p_over = ice->rho * standard_gravity * thk;  // the pressure of the weight of the ice

  p_eff  = p_over - getBasalWaterPressure(thk, bwat, bmr, frac, hmelt_max);
</code>
 */
PetscScalar PISMDefaultYieldStress::basal_water_pressure(PetscScalar thk, PetscScalar bwat,
                                                  PetscScalar bmr, PetscScalar frac,
                                                  PetscScalar hmelt_max) {
  if (bwat > hmelt_max + 1.0e-6) {
    PetscPrintf(grid.com,
                "PISM ERROR:  bwat = %12.8f exceeds hmelt_max = %12.8f\n"
                "  in IceModel::getBasalWaterPressure()\n",bwat,hmelt_max);
    PISMEnd();
  }

  // the model; note 0 <= p_pw <= frac * p_overburden
  // because  0 <= bwat <= hmelt_max
  const PetscScalar p_overburden = p.ice_density * p.standard_gravity * thk; // FIXME task #7297
  PetscScalar p_pw = frac * (bwat / hmelt_max) * p_overburden;

  if (p.usebmr) {
    // add to pressure from instantaneous basal melt rate;
    //   note  (additional) <= (1.0 - frac) * p_overburden so
    //   0 <= p_pw <= p_overburden
    p_pw += ( 1.0 - exp( - PetscMax(0.0,bmr) / p.bmr_scale ) )
      * (1.0 - frac) * p_overburden;
  }

  if (p.usethkeff) {
    // ice thickness is surrogate for distance to margin; near margin the till
    //   is presumably better drained so we reduce the water pressure
    if (thk < p.thkeff_H_high) {
      if (thk <= p.thkeff_H_low) {
        p_pw *= p.thkeff_reduce;
      } else {
        // case Hlow < thk < Hhigh; use linear to connect (Hlow, reduced * p_pw)
        //   to (Hhigh, 1.0 * p_w)
        p_pw *= p.thkeff_reduce
          + (1.0 - p.thkeff_reduce)
          * (thk - p.thkeff_H_low) / (p.thkeff_H_high - p.thkeff_H_low);
      }
    }
  }

  return p_pw;
}

PYS_bwp::PYS_bwp(PISMDefaultYieldStress *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMDefaultYieldStress>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("bwp", grid);
  set_attrs("subglacial (pore) water pressure", "", "Pa", "Pa", 0);
  vars[0].set("_FillValue", -0.01);
  vars[0].set("valid_min", 0);
}

/*!
  \f[p_w = \alpha\, \frac{w}{w_{\text{max}}}\, \rho\, g\, H,\f]
  where 

  - \f$\alpha\f$ is the till pore water fraction (till_pw_fraction),
  - \f$w\f$ is the effective thickness of subglacial melt water (bwat)
  - \f$w_{\text{max}}\f$ is the maximum allowed value for \f$w\f$ (hmelt_max),
  - \f$\rho\f$ is the ice density (ice_density)
  - \f$H\f$ is the ice thickness (thk)

Result is set to invalid (_FillValue) where the ice is floating, there being
no meaning to the above calculation.
 */
PetscErrorCode PYS_bwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bwp", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  const PetscScalar
    alpha     = model->config.get("till_pw_fraction"),
    wmax      = model->config.get("hmelt_max"),
    fillval   = -0.01;

  ierr = model->ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = model->basal_water_thickness->begin_access(); CHKERRQ(ierr);
  ierr = model->basal_melt_rate->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal thk = (*model->ice_thickness)(i, j);
      if (thk > 0.0) {
        (*result)(i,j) = model->basal_water_pressure(thk, // FIXME task #7297
                                                     (*model->basal_water_thickness)(i,j),
                                                     (*model->basal_melt_rate)(i,j),
                                                     alpha, wmax);
      } else { // put negative value below valid range
        (*result)(i,j) = fillval;
      }
    }
  }
  ierr = model->ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = model->basal_water_thickness->end_access(); CHKERRQ(ierr);
  ierr = model->basal_melt_rate->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  MaskQuery m(*model->mask);

  ierr = m.fill_where_floating(*result, fillval); CHKERRQ(ierr);
  
  output = result;
  return 0;
}
