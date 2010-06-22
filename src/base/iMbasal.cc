// Copyright (C) 2004--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"

/*
This file collects procedures related to SSA-as-sliding law in grounded
areas.  IceModel::basalVelocitySIA() is in iMsia.cc (and is not recommended,
generally).
*/

/*** for ice stream regions (MASK_DRAGGING): ***/
PetscScalar IceModel::basalDragx(PetscScalar **tauc,
                                 PISMVector2 **uv,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], uv[i][j].u, uv[i][j].v);
}

PetscScalar IceModel::basalDragy(PetscScalar **tauc,
                                 PISMVector2 **uv,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], uv[i][j].u, uv[i][j].v);
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
relation (\ref Paterson, \ref SchoofStream),
    \f[   \tau_c = c_0 + (\tan \varphi) N, \f]
where \f$N = \rho g H - p_w\f$ is the effective pressure on the till (see
getBasalWaterPressure()), the determination of the till friction angle
\f$\varphi(x,y)\f$  is important.  The till friction angle is assumed to be a
time-independent factor which describes the strength of the unsaturated "dry"
till material.  Thus it is assumed to change more slowly than the basal water 
pressure.  It follows that it changes more slowly than the yield stress and/or
the basal shear stress.

The current method determines the map of till friction angle \c vtillphi
according to options.  Option \c -surf_vel_to_phi causes a call to
invertSurfaceVelocities(), <i>but this method is a stub, and currently
unimplemented</i>.  Option \c -topg_to_phi causes call to
computePhiFromBedElevation().  If neither option is given, the current method
leaves \c vtillphi unchanged, and thus either in its read-in-from-file state or 
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
\f$q=0\f$, which is the purely plastic case; in that case there is \e no direct
relation between the yield stress and the sliding velocity, and the difference 
between the driving stress and the yield stress is entirely held by the membrane
stresses.  (There is also no singular mathematical operation as \f$A^q = A^0 = 1\f$.)
 */
PetscErrorCode IceModel::initBasalTillModel() {
  PetscErrorCode ierr;

  PetscScalar pseudo_plastic_q = config.get("pseudo_plastic_q"),
    pseudo_plastic_uthreshold = config.get("pseudo_plastic_uthreshold") / secpera,
    plastic_regularization = config.get("plastic_regularization") / secpera;

  bool do_pseudo_plastic_till = config.get_flag("do_pseudo_plastic_till"),
    use_ssa_velocity = config.get_flag("use_ssa_velocity");

  if (basal == NULL)
    basal = new IceBasalResistancePlasticLaw(plastic_regularization, do_pseudo_plastic_till, 
                                 pseudo_plastic_q, pseudo_plastic_uthreshold);
  
  if (use_ssa_velocity) {
    ierr = basal->printInfo(3,grid.com); CHKERRQ(ierr);
  }

  ierr = vtauc.set(config.get("default_tauc")); CHKERRQ(ierr);

  bool topgphiSet,svphiSet;
  string filename;
  ierr = PetscOptionsBegin(grid.com, "", "Options controlling the basal till model", ""); CHKERRQ(ierr);
  {
    // initialize till friction angle (vtillphi) from options
    ierr = PISMOptionsIsSet("-topg_to_phi", 
      "Use the till friction angle parameterization", topgphiSet); CHKERRQ(ierr);
    ierr = PISMOptionsString("-surf_vel_to_phi", 
      "Specifies the file containing surface velocities to invert",
			filename, svphiSet); CHKERRQ(ierr);
			
	  bool scaleSet;
	  double slidescale;
    ierr = PISMOptionsReal("-sliding_scale", 
      "Divides pseudo-plastic tauc (yield stress) by given factor; this would increase sliding by given factor in absence of membrane stresses",
      slidescale,scaleSet); CHKERRQ(ierr);
    if (scaleSet) { // only modify config if option set; otherwise leave alone
      if (slidescale > 0.0) {
        ierr = verbPrintf(2, grid.com, 
          "option -sliding_scale read; pseudo yield stress tauc will be divided by %.4f to\n"
          "  cause notional sliding speed-up by given factor %.4f ...\n",
          pow(slidescale,pseudo_plastic_q),slidescale);  CHKERRQ(ierr);
        config.set("sliding_scale_factor_reduces_tauc",slidescale);
      } else {
        ierr = verbPrintf(1, grid.com, 
          "PISM WARNING: negative or zero value given for option -sliding_scale ignored\n");
          CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (svphiSet) {
    SETERRQ(1,"options -surf_vel_to_phi is not available ... ENDING ...\n");
  }

  if (topgphiSet) {
    ierr = verbPrintf(2, grid.com, 
      "option -topg_to_phi seen; creating till friction angle map from bed elev ...\n");
      CHKERRQ(ierr);
    // note option -topg_to_phi will be read again to get comma separated array of parameters
    ierr = computePhiFromBedElevation(); CHKERRQ(ierr);
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
PetscErrorCode IceModel::computePhiFromBedElevation() {

  PetscErrorCode ierr;

  PetscInt    Nparam=5;
  PetscReal   inarray[5] = {5.0, 15.0, -1000.0, 1000.0, 10.0};

  // read comma-separated array of zero to five values
  PetscTruth  topgphiSet;
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-topg_to_phi", inarray, &Nparam, &topgphiSet);
     CHKERRQ(ierr);
  if (topgphiSet != PETSC_TRUE) {
    SETERRQ(1,"HOW DID I GET HERE? ... ending...\n");
  }
  if ((Nparam > 5) || (Nparam < 4)) {
    ierr = verbPrintf(1, grid.com, 
      "PISM ERROR: option -topg_to_phi provided with more than 5 or fewer than 4\n"
      "            arguments ... ENDING ...\n");
      CHKERRQ(ierr);
    PetscEnd();
  }
  PetscReal   phi_min = inarray[0],
              phi_max = inarray[1],
              topg_min = inarray[2],
              topg_max = inarray[3],
              phi_ocean = inarray[4];

  ierr = verbPrintf(2, grid.com, 
      "  till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
      "            /  %5.2f                                 for   topg < %.f\n"
      "      phi = |  %5.2f + (topg - %.f) * (%.2f / %.f)   for   %.f < topg < %.f\n"
      "            \\  %5.2f                                 for   %.f < topg\n",
      phi_min, topg_min,
      phi_min, topg_min, phi_max-phi_min, topg_max - topg_min, topg_min, topg_max,
      phi_max, topg_max);
      CHKERRQ(ierr);
  if (Nparam == 5) {
    ierr = verbPrintf(2, grid.com, 
      "      (also using phi = %5.2f in floating ice or ice free ocean)\n",
      phi_ocean); CHKERRQ(ierr);
  }

  PetscReal slope = (phi_max - phi_min) / (topg_max - topg_min);
  PetscScalar **tillphi, **bed;
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vbed.get_array(bed); CHKERRQ(ierr);
  ierr = vtillphi.get_array(tillphi); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((!vMask.is_floating(i,j)) || (Nparam < 5)) {
        if (bed[i][j] <= topg_min) {
          tillphi[i][j] = phi_min;
        } else if (bed[i][j] >= topg_max) {
          tillphi[i][j] = phi_max;
        } else {
          tillphi[i][j] = phi_min + (bed[i][j] - topg_min) * slope;
        }
      } else {
        tillphi[i][j] = phi_ocean;
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vtillphi.end_access(); CHKERRQ(ierr);

  return 0;
}


//! Compute modeled pressure in subglacial liquid water using thickness of water layer, and possibly the melt rate, at the base.
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
PetscScalar IceModel::getBasalWaterPressure(PetscScalar thk, PetscScalar bwat,
				            PetscScalar bmr, PetscScalar frac,
				            PetscScalar hmelt_max) const {

  if (bwat > hmelt_max) {
    PetscPrintf(grid.com,
      "PISM ERROR:  bwat exceeds hmelt_max in IceModel::getBasalWaterPressure()\n");
    PetscEnd();
  }

  static const bool
    usebmr        = config.get_flag("bmr_enhance_basal_water_pressure"),
    usethkeff     = config.get_flag("thk_eff_basal_water_pressure");
  static const PetscScalar
    bmr_scale     = config.get("bmr_enhance_scale"),
    thkeff_reduce = config.get("thk_eff_reduced"),
    thkeff_H_high = config.get("thk_eff_H_high"),
    thkeff_H_low  = config.get("thk_eff_H_low");

  // the model; note  0 <= p_pw <= frac * p_overburden
  //   because  0 <= bwat <= hmelt_max
  const PetscScalar p_overburden = ice->rho * standard_gravity * thk;
  PetscScalar  p_pw = frac * (bwat / hmelt_max) * p_overburden;

  if (usebmr) {
    // add to pressure from instantaneous basal melt rate;
    //   note  (additional) <= (1.0 - frac) * p_overburden so
    //   0 <= p_pw <= p_overburden
    p_pw += ( 1.0 - exp( - PetscMax(0.0,bmr) / bmr_scale ) )
            * (1.0 - frac) * p_overburden;
  }

  if (usethkeff) {
    // ice thickness is surrogate for distance to margin; near margin the till
    //   is presumably better drained so we reduce the water pressure
    if (thk < thkeff_H_high) {
      if (thk <= thkeff_H_low) {
        p_pw *= thkeff_reduce;
      } else {
        // case Hlow < thk < Hhigh; use linear to connect (Hlow, reduced * p_pw)
        //   to (Hhigh, 1.0 * p_w)
        p_pw *= thkeff_reduce
                + (1.0 - thkeff_reduce)
                    * (thk - thkeff_H_low) / (thkeff_H_high - thkeff_H_low);
      }
    }
  }

  return p_pw;
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
PetscErrorCode IceModel::updateYieldStressUsingBasalWater() {
  PetscErrorCode  ierr;

  bool use_ssa_when_grounded = config.get_flag("use_ssa_when_grounded");
  // only makes sense when use_ssa_when_grounded == TRUE
  if (use_ssa_when_grounded == PETSC_FALSE) {
    SETERRQ(1,"use_ssa_when_grounded == PETSC_FALSE but updateYieldStressFromHmelt() called");
  }

  if (holdTillYieldStress == PETSC_FALSE) { // usual case: use Hmelt to determine tauc
    PetscScalar till_pw_fraction = config.get("till_pw_fraction"),
      till_c_0 = config.get("till_c_0") * 1e3, // convert from kPa to Pa
      till_mu = tan((pi/180.0)*config.get("default_till_phi")),
      hmelt_max = config.get("hmelt_max");

    ierr =    vMask.begin_access(); CHKERRQ(ierr);
    ierr =    vtauc.begin_access(); CHKERRQ(ierr);
    ierr =       vH.begin_access(); CHKERRQ(ierr);
    ierr =   vHmelt.begin_access(); CHKERRQ(ierr);
    ierr =     vbmr.begin_access(); CHKERRQ(ierr);
    ierr = vtillphi.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (vMask.is_floating(i,j)) {
          vtauc(i,j) = 0.0;  
        } else if (vH(i,j) == 0.0) {
          vtauc(i,j) = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar if no ice
        } else { // grounded and there is some ice
          const PetscScalar
            p_over = ice->rho * standard_gravity * vH(i,j),
            p_w    = getBasalWaterPressure(vH(i,j), vHmelt(i,j),
                       vbmr(i,j), till_pw_fraction, hmelt_max),
            N      = p_over - p_w;  // effective pressure on till
          if (useConstantTillPhi == PETSC_TRUE) {
            vtauc(i,j) = till_c_0 + N * till_mu;
          } else {
            vtauc(i,j) = till_c_0 + N * tan((pi/180.0) * vtillphi(i,j));
          }
        }
      }
    }
    ierr =    vMask.end_access(); CHKERRQ(ierr);
    ierr =    vtauc.end_access(); CHKERRQ(ierr);
    ierr =       vH.end_access(); CHKERRQ(ierr);
    ierr = vtillphi.end_access(); CHKERRQ(ierr);
    ierr =     vbmr.end_access(); CHKERRQ(ierr);
    ierr =   vHmelt.end_access(); CHKERRQ(ierr);
  }

  // scale tauc if desired
  double A = config.get("sliding_scale_factor_reduces_tauc");
  if (A > 0.0) {
    const PetscScalar q = config.get("pseudo_plastic_q");
    vtauc.scale(1.0 / pow(A,q));
  }

  return 0;
}



//! Apply explicit time step for pure diffusion to basal layer of melt water.
/*!
See preprint \ref BBssasliding .

Uses vWork2d[0] to temporarily store new values for Hmelt.
 */
PetscErrorCode IceModel::diffuseHmelt() {
  PetscErrorCode  ierr;

  const PetscScalar
    L = config.get("hmelt_diffusion_distance"),
    diffusion_time = config.get("hmelt_diffusion_time") * secpera; // convert to seconds

  // diffusion constant K in u_t = K \nabla^2 u is chosen so that fundmental
  //   solution has standard deviation \sigma = 20 km at time t = 1000 yrs;
  //   2 \sigma^2 = 4 K t
  const PetscScalar K = L * L / (2.0 * diffusion_time),
                    Rx = K * dtTempAge / (grid.dx * grid.dx),
                    Ry = K * dtTempAge / (grid.dy * grid.dy);

  // NOTE: restriction that
  //    1 - 2 R_x - 2 R_y \ge 0
  // is a maximum principle restriction; therefore new Hmelt will be between
  // zero and hmelt_max if old Hmelt has that property
  const PetscScalar oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  if (oneM4R <= 0.0) {
    SETERRQ(1,
       "diffuseHmelt() has 1 - 2Rx - 2Ry <= 0 so explicit method for diffusion unstable\n"
       "  (timestep restriction believed so rare that is not part of adaptive scheme)");
  }

  // communicate ghosted values so neighbors are valid (temperatureStep and
  // enthalpyAndDrainageStep modify vHmelt, but do not update ghosts).
  ierr = vHmelt.beginGhostComm(); CHKERRQ(ierr);
  ierr = vHmelt.endGhostComm(); CHKERRQ(ierr);

  PetscScalar **Hmelt, **Hmeltnew; 
  ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(Hmeltnew); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Hmeltnew[i][j] = oneM4R * Hmelt[i][j]
                       + Rx * (Hmelt[i+1][j] + Hmelt[i-1][j])
                       + Ry * (Hmelt[i][j+1] + Hmelt[i][j-1]);
    }
  }
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  // finally copy new into vHmelt (and communicate ghosted values at the same time)
  ierr = vWork2d[0].beginGhostComm(vHmelt); CHKERRQ(ierr);
  ierr = vWork2d[0].endGhostComm(vHmelt); CHKERRQ(ierr);

  return 0;
}

