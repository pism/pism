// Copyright (C) 2008--2009 Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include "nc_util.hh"
#include "LocalInterpCtx.hh"
#include "iceModelVec.hh"
#include "iceModel.hh"

/*  for example usage see comments in examples/pst/convert_to_inv.sh; note 
P0A.nc is at ftp://ftp.gi.alaska.edu/pub/bueler/P0A.nc  */

/*
INVERSE MODEL; with and w/o regularization:

$ pismr -ssa -super -plastic -i inv_me.nc -y 1 -pseudo_plastic_q 0.25 -surf_vel_to_phi inv_me.nc \
   -inv_write_fields foo.nc -o inv_result.nc

$ pismr -ssa -super -plastic -i inv_me.nc -y 1 -pseudo_plastic_q 0.25 -surf_vel_to_phi inv_me.nc \
   -inv_write_fields foo_noreg.nc -inv_reg_eps 0.0 -o inv_result_noreg.nc

INVERSE MODEL options:

   -inv_reg_eps              argument of 0.0 turns off regularization; default = 1.0e6
   -inv_phi_min              default = 5.0
   -inv_phi_max              default = 15.0
   -inv_write_fields foo.nc  write several fields associated to inverse model to foo.nc
*/


//! Invert given ice surface velocities to find till friction angle using a pseudo-plastic model.
/*!
Expects to find variables \c uvelsurf and \c vvelsurf, namely x and y components
of observed surface velocity, in NetCDF file \c filename.  Calls 
readObservedSurfVels() for this purpose.  If found, reads them.

Enough space for ten IceModelVec2's is allocated.  In particular,
this procedure creates and destroys the InverseModelCtx (instance IceModel::inv) members.

The first goal is to construct a mask showing where velocities are present; this is
also done by readObservedSurfVels().  A point is marked as not having a valid
observed velocity (=0) if either x or y component has value outside the valid range
specified in the NetCDF file for that variable.  (Thus the creator of the input file
is in charge of setting \c valid_range or \c valid_min and \c valid_max attributes
for both variables \c uvelsurf and \c vvelsurf.  \c _FillValue, if set, should
be outside the valid range.)  A valid point is marked as 2 if there is stencil width for
differentiating the observed velocity, and 1 otherwise.

Then calls computeSIASurfaceVelocity() to get the value of the surface
velocity which would be computed from SIA flow model for the current
geometry and temperature.

If <tt>-super</tt> is also set, we call the method computeFofVforInverse()
which computes the factor \f$f(|\mathbf{v}|)\f$ used in combining the SIA 
and SSA velocities \lo\cite[equation (21)]{BBssasliding}\elo.  If
<tt>-super</tt> is not set, we do not call computeFofVforInverse(), and
instead just set \f$f(|\mathbf{v}|)=0\f$, so the observed surface velocity
is treated as all described by the SSA model.

This method then calls three steps in turn:

- removeSIApart() to compute the basal sliding velocity \f$\mathbf{u}_b\f$ by removing the SIA
modeled shear from the surface velocity,

- computeBasalShearFromSSA() to compute the basal shear stress \f$\tau_b\f$ which solves 
the SSA stress balance in which the basal sliding velocity \f$\mathbf{u}_b\f$ is used as the
depth independent velocity,

- \e either computeTFAFromBasalShear() to compute the till friction angle \f$\phi\f$
for the given basal shear stress \f$\tau_b\f$ and basal sliding velocity 
\f$\mathbf{u}_b\f$ using a regularized objective, \e or 
computeTFAFromBasalShearNoReg() with no regularization if the option
<tt>-inv_reg_eps 0.0</tt> is used to turn off the regularization.

The result is a field of values for till friction angle in a pseudo-plastic till model.  
PISM can do the forward time-stepping model using the resulting till friction angle
field.

This model includes purely-plastic and linearly-viscous cases, but is expected to 
be most effective in intermediate pseudo-plastic cases.  There is no 
smoothing in this initial implementation.  This method must be called after 
initialization is essentially completed.  It uses geometry, temperature field, 
and map of effective thickness of basal till water in the inversion.

This method also reads user option <tt>-inv_write_fields bar.nc</tt>.  If
this option is found, writes most intermediate fields from the inverse model
computation to bar.nc.
 */
PetscErrorCode IceModel::invertSurfaceVelocities(const char *filename) {
  PetscErrorCode ierr;
  
  // read options
  PetscTruth  invfieldsSet;
  PetscScalar invPhiMax = 15.0, 
              invPhiMin = 5.0,
              invRegEps = 1.0e23;
  char invfieldsname[PETSC_MAX_PATH_LEN];

  if (doPseudoPlasticTill == PETSC_FALSE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: inverse model and invertSurfaceVelocites() should only\n"
       "  be called with q > 0.0 in pseudo-plastic model;  here is \n"
       "  PlasticBasalType::printInfo() output:\n");
       CHKERRQ(ierr);
    ierr = basal->printInfo(1,grid.com); CHKERRQ(ierr);
    ierr = verbPrintf(1, grid.com, "  CONTINUING.  May crash!!\n"); CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-inv_phi_min", &invPhiMin, PETSC_NULL);
           CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-inv_phi_max", &invPhiMax, PETSC_NULL);
           CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,   "-inv_reg_eps", &invRegEps, PETSC_NULL);
           CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-inv_write_fields", invfieldsname, 
           PETSC_MAX_PATH_LEN, &invfieldsSet); CHKERRQ(ierr);

  // allocate
  ierr = createInvFields(); CHKERRQ(ierr);

  // if user will want fields written, prepare the file now
  if (invfieldsSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
             "  preparing file %s to write inverse computation fields ...\n",
             invfieldsname); CHKERRQ(ierr);
    NCTool nc(&grid);
    ierr = nc.open_for_writing(invfieldsname,PETSC_TRUE); CHKERRQ(ierr);
    ierr = nc.append_time(grid.year * secpera); CHKERRQ(ierr);
    ierr = nc.write_history("option -inv_write_fields read"); CHKERRQ(ierr);
    ierr = nc.write_polar_stereographic(psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr);
    ierr = nc.write_global_attrs(PETSC_FALSE, "CF-1.3"); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
  } else {
    strcpy(invfieldsname,""); // make sure empty
  }

  // copy till phi to old
  ierr = inv.oldtillphi->copy_from(vtillphi); CHKERRQ(ierr);
  
  // read in surface velocity
  ierr = verbPrintf(2, grid.com, 
     "  reading observed surface velocities from file %s ...\n"
     "  computing mask for where observed surfaced velocities are present ...\n",
     filename); CHKERRQ(ierr);
  ierr = readObservedSurfVels(filename); CHKERRQ(ierr);

#if 1
  // smooth input velocities (FIXME: keeping mask as is, for now)
  ierr = verbPrintf(2, grid.com,
     "  smoothing observed velocities ...\n"); CHKERRQ(ierr);
  const PetscInt smoothingPasses = 3;
  ierr = smoothObservedSurfVels(smoothingPasses); CHKERRQ(ierr);
#endif

  // compute the surface velocity which the nonsliding SIA predicts
  ierr = verbPrintf(2, grid.com, 
    "  computing SIA surface velocity (for removal from observed before inverting w SSA ...\n");
    CHKERRQ(ierr);
  ierr = computeSIASurfaceVelocity(); CHKERRQ(ierr);

  // compute f(|v|) factor, or set to constant if -super not used
  if (doSuperpose == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
       "  flag doSuperpose (option -super) seen;  computing f(|v|) from observed\n"
       "    surface velocities by solving transcendental equations ...\n"); CHKERRQ(ierr);
    ierr = computeFofVforInverse(); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
       "  flag doSuperpose (option -super) NOT seen;  setting f(|v|) to 0.0, so NONE\n"
       "    of SIA velocity is removed from observed surface velocity ...\n");
       CHKERRQ(ierr);
    ierr = inv.fofv->set(0.0); CHKERRQ(ierr);
  }

  // major step of inverse model: remove SIA
  ierr = verbPrintf(2, grid.com, 
           "  removing SIA-computed (vertical shear) part of surface velocity ...\n");
           CHKERRQ(ierr);
  ierr = removeSIApart();    CHKERRQ(ierr);

  // major step of inverse model: apply SSA operator
  ierr = verbPrintf(2, grid.com, 
           "  computing basal shear stress: applying SSA differential operator to velocity ...\n");
           CHKERRQ(ierr);
  ierr = computeBasalShearFromSSA();   CHKERRQ(ierr);

  // major step of inverse model: determine till friction angle from basal shear
  ierr = verbPrintf(2, grid.com, 
           "  computing till friction angle using (pseudo-)plastic model and Mohr-Coulomb criterion ...\n"
           "    getting effective pressure on till ...\n"); 
           CHKERRQ(ierr);
  ierr = getEffectivePressureForInverse(); CHKERRQ(ierr); // need N regardless of regularization
  if (invRegEps == 0.0) {
    ierr = verbPrintf(2, grid.com, 
           "    NOT regularizing mu = tan(phi) computation (epsilon = 0.0) ...\n");
           CHKERRQ(ierr);
    ierr = computeTFAFromBasalShearNoReg(invPhiMin,invPhiMax); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
           "    regularizing mu = tan(phi) computation using epsilon = %.3e ...\n",
           invRegEps); CHKERRQ(ierr);
    ierr = computeTFAFromBasalShear(invPhiMin,invPhiMax,invRegEps,invfieldsname); CHKERRQ(ierr);
  }

  // write out stored inverse info; mostly for debug
  if (invfieldsSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
             "  writing various fields from inverse model computation to file %s ...\n",
             invfieldsname); CHKERRQ(ierr);
    ierr = writeInvFields(invfieldsname); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com, 
     "done with ad hoc inverse model; running forward with new till friction angle map ...\n",
     filename); CHKERRQ(ierr);
  
  ierr = destroyInvFields(); CHKERRQ(ierr);
  return 0;
}


//! Allocate IceModelVec2's for inverse model and set metadata.
PetscErrorCode IceModel::createInvFields() {
  PetscErrorCode ierr;
  inv.usIn = new IceModelVec2;
  ierr = inv.usIn->create(grid, "uvelsurf", true); CHKERRQ(ierr);
  ierr = inv.usIn->set_attrs(
     "inverse_input", 
     "x component of observed velocity of ice at ice surface","m s-1", NULL); CHKERRQ(ierr);
  // only needed for std out range message at read:
  ierr = inv.usIn->set_glaciological_units("m year-1"); CHKERRQ(ierr);

  inv.vsIn = new IceModelVec2;
  ierr = inv.vsIn->create(grid, "vvelsurf", true); CHKERRQ(ierr);
  ierr = inv.vsIn->set_attrs(
     "inverse_input", 
     "y component of observed velocity of ice at ice surface","m s-1", NULL); CHKERRQ(ierr);
  // only needed for display message at read:
  ierr = inv.vsIn->set_glaciological_units("m year-1"); CHKERRQ(ierr);

  inv.invMask = new IceModelVec2;
  ierr = inv.invMask->create(grid, "invMask", true); CHKERRQ(ierr);
  ierr = inv.invMask->set_attrs(
     "inverse_output", 
     "mask for presence of observed velocity at ice surface", NULL, NULL); CHKERRQ(ierr);

  inv.usSIA = new IceModelVec2;
  ierr = inv.usSIA->create(grid, "uvelsurfSIA", false); CHKERRQ(ierr);  // global
  ierr = inv.usSIA->set_attrs(
     "inverse_output", 
     "x component of ice surface velocity predicted by non-sliding SIA",
     "m s-1", NULL); CHKERRQ(ierr);

  inv.vsSIA = new IceModelVec2;
  ierr = inv.vsSIA->create(grid, "vvelsurfSIA", false); CHKERRQ(ierr);  // global
  ierr = inv.vsSIA->set_attrs(
     "inverse_output",
     "y component of ice surface velocity predicted by non-sliding SIA",
     "m s-1", NULL); CHKERRQ(ierr);

  inv.taubxComputed = new IceModelVec2;
  ierr = inv.taubxComputed->create(grid, "taubxOUT", false);  // global
  ierr = inv.taubxComputed->set_attrs(
     "inverse_output", 
     "inverse-model-computed x component of basal shear stress", 
     "Pa", NULL); CHKERRQ(ierr);

  inv.taubyComputed = new IceModelVec2;
  ierr = inv.taubyComputed->create(grid, "taubyOUT", false);  // global
  ierr = inv.taubyComputed->set_attrs(
     "inverse_output",
     "inverse-model-computed y component of basal shear stress", 
     "Pa", NULL); CHKERRQ(ierr);

  inv.fofv = new IceModelVec2;
  ierr = inv.fofv->create(grid, "fofv", false);  // global
  ierr = inv.fofv->set_attrs(
     "inverse_output",
     "inverse-model-computed fraction of velocity from SIA in hybrid model", 
     NULL, NULL); CHKERRQ(ierr);

  inv.oldtillphi = new IceModelVec2;
  ierr = inv.oldtillphi->create(grid, "tillphiOLD", true);  // local like vtillphi
  ierr = inv.oldtillphi->set_attrs(
     "inverse_output",
     "till friction angle at start of inverse model", 
     NULL, NULL); CHKERRQ(ierr);

  inv.effPressureN = new IceModelVec2;
  ierr = inv.effPressureN->create(grid, "effpress", false);
  ierr = inv.effPressureN->set_attrs(
     "inverse_output",
     "effective pressure on till material", 
     "Pa", NULL); CHKERRQ(ierr);

  return 0;
}


//! De-allocate IceModelVec2's for inverse model.
PetscErrorCode IceModel::destroyInvFields() {
  PetscErrorCode ierr;
  ierr = inv.usIn->destroy(); CHKERRQ(ierr);
  delete inv.usIn;
  inv.usIn = PETSC_NULL;
  ierr = inv.vsIn->destroy(); CHKERRQ(ierr);
  delete inv.vsIn;
  inv.vsIn = PETSC_NULL;
  ierr = inv.invMask->destroy(); CHKERRQ(ierr);
  delete inv.invMask;
  inv.invMask = PETSC_NULL;
  ierr = inv.usSIA->destroy(); CHKERRQ(ierr);
  delete inv.usSIA;
  inv.usSIA = PETSC_NULL;
  ierr = inv.vsSIA->destroy(); CHKERRQ(ierr);
  delete inv.vsSIA;
  inv.vsSIA = PETSC_NULL;
  ierr = inv.taubxComputed->destroy(); CHKERRQ(ierr);
  delete inv.taubxComputed;
  inv.taubxComputed = PETSC_NULL;
  ierr = inv.taubyComputed->destroy(); CHKERRQ(ierr);
  delete inv.taubyComputed;
  inv.taubyComputed = PETSC_NULL;
  ierr = inv.fofv->destroy(); CHKERRQ(ierr);
  delete inv.fofv;
  inv.fofv = PETSC_NULL;
  ierr = inv.oldtillphi->destroy(); CHKERRQ(ierr);
  delete inv.oldtillphi;
  inv.oldtillphi = PETSC_NULL;
  delete inv.effPressureN;
  inv.effPressureN = PETSC_NULL;
  return 0;
}


//! Read observed surface velocities and build mask for validity.
/*!
Mask is 0 where no observed velocity is present, and 1 or 2 where it is present.
Mask is 2 when a 9 point stencil of valid velocities surrounds the point, so that
the SSA differential operator can be applied.

Validity is determined entirely by querying the input file.  The input file
should have x,y components of velocity called "uvelsurf", "vvelsurf".  These
variables should each have either valid_range or both valid_min and valid_max.
The mask is set to zero where the value is outside of this range.  Presumably
the invalid points have a value which matches _FillValue for that variable, but
all that is required is that the fill value is outside the valid range.
 */
PetscErrorCode IceModel::readObservedSurfVels(const char *filename) {
  PetscErrorCode ierr;
  bool file_exists = false; // check whether file exists
  NCTool nc(&grid);
  grid_info gi;

  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com,"PISM ERROR: Can't open file '%s'.\n",
    filename); CHKERRQ(ierr);    PetscEnd();
  }
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  LocalInterpCtx lic(gi, NULL, NULL, grid); // 2D only
  // read in by regridding; checks for existence of file and variable in file,
  //   and reads valid_ attributes
  ierr = inv.usIn->regrid(filename, lic, true); CHKERRQ(ierr);// it *is* critical
  ierr = inv.vsIn->regrid(filename, lic, true); CHKERRQ(ierr);

  PetscScalar **us, **vs, **imask;

  ierr = inv.usIn->get_array(us);      CHKERRQ(ierr);
  ierr = inv.vsIn->get_array(vs);      CHKERRQ(ierr);
  ierr = inv.invMask->get_array(imask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (inv.usIn->is_valid(us[i][j]) && inv.vsIn->is_valid(vs[i][j])) {
        imask[i][j] = 1.0;
      } else {
        imask[i][j] = 0.0;
      }
    }
  }
  ierr = inv.usIn->end_access(); CHKERRQ(ierr);
  ierr = inv.vsIn->end_access(); CHKERRQ(ierr);
  ierr = inv.invMask->end_access(); CHKERRQ(ierr);

  // need stencil width on all:
  ierr = inv.usIn->beginGhostComm(); CHKERRQ(ierr);
  ierr = inv.vsIn->beginGhostComm(); CHKERRQ(ierr);
  ierr = inv.invMask->beginGhostComm(); CHKERRQ(ierr);
  ierr = inv.usIn->endGhostComm(); CHKERRQ(ierr);
  ierr = inv.vsIn->endGhostComm(); CHKERRQ(ierr);
  ierr = inv.invMask->endGhostComm(); CHKERRQ(ierr);

  ierr = inv.invMask->get_array(imask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // invMask is 2 if SSA differential operator can differentiate 9 point stencil
      if (    (imask[i-1][j-1] > 0.5) && (imask[i][j-1] > 0.5) && (imask[i+1][j-1] > 0.5)
           && (imask[i-1][j] > 0.5)   && (imask[i][j] > 0.5)   && (imask[i+1][j] > 0.5)
           && (imask[i-1][j+1] > 0.5) && (imask[i][j+1] > 0.5) && (imask[i+1][j+1] > 0.5)  ) {
        imask[i][j] = 2.0;
      }
    }
  }
  ierr = inv.invMask->end_access(); CHKERRQ(ierr);

  // again, for the mask
  ierr = inv.invMask->beginGhostComm(); CHKERRQ(ierr);
  ierr = inv.invMask->endGhostComm(); CHKERRQ(ierr);

  return 0;
}


//! Write fields associated to inverse model to a prepared NetCDF file.
PetscErrorCode IceModel::writeInvFields(const char *filename) {
  PetscErrorCode ierr;
  PetscScalar fill_ma  = 2.0 * secpera; // 2.0 m/s is the fill value; huge

  ierr = inv.usIn->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  inv.usIn->write_in_glaciological_units = true;
  ierr = inv.usIn->write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = inv.usIn->write_scalar_attr(filename, "_FillValue", NC_FLOAT, 1, &fill_ma); CHKERRQ(ierr);

  ierr = inv.vsIn->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  inv.vsIn->write_in_glaciological_units = true;
  ierr = inv.vsIn->write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = inv.vsIn->write_scalar_attr(filename, "_FillValue", NC_FLOAT, 1, &fill_ma); CHKERRQ(ierr);

  ierr = getMagnitudeOf2dVectorField(*(inv.usIn), *(inv.usIn), vWork2d[0]); CHKERRQ(ierr);
  ierr = vWork2d[0].set_name("magvelsurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("inverse_output",
             "magnitude of observed velocity of ice at ice surface",
	     "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = inv.invMask->write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = inv.usSIA->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  inv.usSIA->write_in_glaciological_units = true;
  ierr = inv.usSIA->write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = inv.vsSIA->set_glaciological_units("m year-1"); CHKERRQ(ierr);
  inv.vsSIA->write_in_glaciological_units = true;
  ierr = inv.vsSIA->write(filename, NC_FLOAT); CHKERRQ(ierr);

  bool oldwritegu = vubarSSA.write_in_glaciological_units;
  ierr = vubarSSA.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vubarSSA.write_in_glaciological_units = true;
  ierr = vubarSSA.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = vubarSSA.write_scalar_attr(filename, "_FillValue", NC_FLOAT, 1, &fill_ma); CHKERRQ(ierr);
  if (!oldwritegu)   vubarSSA.write_in_glaciological_units = false;

  oldwritegu = vvbarSSA.write_in_glaciological_units;
  ierr = vvbarSSA.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vvbarSSA.write_in_glaciological_units = true;
  ierr = vvbarSSA.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = vvbarSSA.write_scalar_attr(filename, "_FillValue", NC_FLOAT, 1, &fill_ma); CHKERRQ(ierr);
  if (!oldwritegu)   vvbarSSA.write_in_glaciological_units = false;

  ierr = inv.taubxComputed->write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = inv.taubyComputed->write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = inv.effPressureN->write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = getMagnitudeOf2dVectorField(*(inv.taubxComputed),*(inv.taubyComputed),
                                     vWork2d[0]); CHKERRQ(ierr);
  ierr = vWork2d[0].set_name("magtaubComputed"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("inverse_output",
             "magnitude of basal shear stress applied at base of ice",
	     "Pa", NULL); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = false;
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = inv.fofv->write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = inv.fofv->write_text_attr(filename, "more_info", 
          "value of 1 means velocity is all SIA, value of 0 means velocity is all SSA");
          CHKERRQ(ierr);

  // input till phi
  ierr = inv.oldtillphi->write(filename, NC_FLOAT); CHKERRQ(ierr);

  // output till phi
  ierr = vtillphi.write(filename, NC_FLOAT); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::smoothObservedSurfVels(const PetscInt passes) {
  PetscErrorCode ierr;
  PetscScalar **us, **vs, **usNEW, **vsNEW, **imask;

  // following is [0.35 0.30 0.35] outer product w itself; ref conversation with Orion L.;
  //   retrieved from src/base/iMvelocity.cc in r174
  const PetscScalar c[3][3] = {{0.1225, 0.105,  0.1225},
                               {0.105,  0.09,   0.105},
                               {0.1225, 0.105,  0.1225}};
  // note communication on inv.usIn, inv.vsIn should have happened in readObservedSurfVels()
  ierr = verbPrintf(2, grid.com, "    smoothing pass: "); CHKERRQ(ierr);
  for (PetscInt m=0; m<passes; ++m) {
    ierr = verbPrintf(2, grid.com, " %d ...", m+1); CHKERRQ(ierr);
    ierr =  inv.usIn->get_array(us);    CHKERRQ(ierr);
    ierr =  inv.vsIn->get_array(vs);    CHKERRQ(ierr);
    ierr = vWork2d[0].get_array(usNEW); CHKERRQ(ierr);
    ierr = vWork2d[1].get_array(vsNEW); CHKERRQ(ierr);
    ierr = inv.invMask->get_array(imask); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (imask[i][j] > 1.5) {
          usNEW[i][j] =  c[0][0]*us[i-1][j-1] + c[0][1]*us[i-1][j] + c[0][2]*us[i-1][j+1]
                       + c[1][0]*us[i][j-1]   + c[1][1]*us[i][j]   + c[1][2]*us[i][j+1]
                       + c[2][0]*us[i+1][j-1] + c[2][1]*us[i+1][j] + c[2][2]*us[i+1][j+1];
          vsNEW[i][j] =  c[0][0]*vs[i-1][j-1] + c[0][1]*vs[i-1][j] + c[0][2]*vs[i-1][j+1]
                       + c[1][0]*vs[i][j-1]   + c[1][1]*vs[i][j]   + c[1][2]*vs[i][j+1]
                       + c[2][0]*vs[i+1][j-1] + c[2][1]*vs[i+1][j] + c[2][2]*vs[i+1][j+1];
        } else { // no change if sans stencil width
          usNEW[i][j] = us[i][j];
          vsNEW[i][j] = vs[i][j];
        }
      }
    }
    ierr =    inv.usIn->end_access(); CHKERRQ(ierr);
    ierr =    inv.vsIn->end_access(); CHKERRQ(ierr);
    ierr =   vWork2d[0].end_access(); CHKERRQ(ierr);
    ierr =   vWork2d[1].end_access(); CHKERRQ(ierr);
    ierr = inv.invMask->end_access(); CHKERRQ(ierr);

    ierr = vWork2d[0].beginGhostComm(*(inv.usIn)); CHKERRQ(ierr);
    ierr = vWork2d[1].beginGhostComm(*(inv.vsIn)); CHKERRQ(ierr);
    ierr = vWork2d[0].endGhostComm(); CHKERRQ(ierr);
    ierr = vWork2d[1].endGhostComm(); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2, grid.com, "done ...\n"); CHKERRQ(ierr);
  return 0;
}

//! Compute the surface velocity predicted by the nonsliding SIA.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

The nonsliding thermomechanically coupled SIA equations are used.

This procedure calls surfaceGradientSIA(), velocitySIAStaggered(), and
horizontalVelocitySIARegular().  Therefore these IceModel members get modified:

- IceModelVec2 vuvbar[2]
- IceModelVec2 vWork2d[0..5]
- IceModelVec2 vub  (set to zero)
- IceModelVec2 vvb  (set to zero)
- IceModelVec3 Istag3[2] 
- IceModelVec3 Sigmastag3[2].
- IceModelVec3 vu
- IceModelVec3 vv
 
Only \c vu and \c vv are directly used here, however.  They are evaluated for 
their surface velocity.
 */
PetscErrorCode IceModel::computeSIASurfaceVelocity() {
  PetscErrorCode ierr;
  PetscScalar    OLDmuSliding = muSliding;
  
  muSliding = 0.0;  // turn off SIA type sliding even if wanted for forward model (a bad idea ...)
  
  ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here ...
  // surface gradient temporarily stored in vWork2d[0..3] ...

  // communicate h_x[o], h_y[o] on staggered for basalSIA() and horizontalVelocitySIARegular()
  for (PetscInt k=0; k<4; ++k) { 
    ierr = vWork2d[k].beginGhostComm(); CHKERRQ(ierr);
    ierr = vWork2d[k].endGhostComm(); CHKERRQ(ierr);
  }

  // compute the vertically integrated velocity from the nonsliding SIA
  //   (namely vuvbar[2], which is (ubar,vbar) on staggered grid)
  ierr = velocitySIAStaggered(); CHKERRQ(ierr);

  // communicate vuvbar[01] for velocities2DSIAToRegular():
  ierr = vuvbar[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[0].endGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[1].endGhostComm(); CHKERRQ(ierr);

  // we have the 2d vertically-averaged velocity on the staggered grid, but
  //   we actually need the surface velocity computed by the nonsliding SIA
  //   on the regular grid, so now we call horizontalVelocitySIARegular();
  // re communication pattern: compare code in velocity()
  ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
  ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);
  ierr = vub.set(0.0); CHKERRQ(ierr);  // no sliding in SIA
  ierr = vvb.set(0.0); CHKERRQ(ierr);
  ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);

  // now grab surface value of 3D velocity fields
  ierr = u3.beginGhostComm(); CHKERRQ(ierr);
  ierr = v3.beginGhostComm(); CHKERRQ(ierr);
  ierr = u3.endGhostComm(); CHKERRQ(ierr);
  ierr = v3.endGhostComm(); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(*(inv.usSIA), vH); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(*(inv.vsSIA), vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  
  muSliding = OLDmuSliding;
  return 0;
}


//! Computes G(x) and its derivative needed by computeFofVforInverse().
/*!
See comments for computeFofVforInverse() for formula for \f$G(x)\f$.
 */
PetscErrorCode IceModel::getGforInverse(const PetscScalar x, 
                   const PetscScalar UsuSIAdiffsqr, const PetscScalar UsuSIAdiffdotuSIA,
                   const PetscScalar uSIAsqr, PetscScalar &G, PetscScalar &Gprime) {

  const PetscScalar  v0sqr  = PetscSqr(100.0 / secpera),
                     outC   = 2.0 / pi,
                     F      = outC * atan(x / v0sqr),
                     Fprime = outC * (1.0 / (1.0 + PetscSqr(x / v0sqr))) / v0sqr;
  G      = (x - uSIAsqr) * F * F - 2.0 * F * UsuSIAdiffdotuSIA - UsuSIAdiffsqr;
  Gprime = F * F + 2.0 * Fprime * ( (x - uSIAsqr) * F - UsuSIAdiffdotuSIA );
  return 0;
}


//! Given observed surface velocity and SIA computed surface velocity, find f(|v|) factor.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

The scalar field which gets computed is called \f$f(|\mathbf{v}|)\f$ in 
\lo\cite[equations (21) and (22)]{BBssasliding}\elo.  Let 
\f$\mathbf{U}_s\f$ be the observed surface velocity.  At each point
in the map plane we seek \f$f(|\mathbf{v}|)\f$ in the equation
  \f[ \mathbf{U}_s = f(|\mathbf{v}|) \mathbf{u}_s
                   + \left(1 - f(|\mathbf{v}|)\right) \mathbf{v} \f]
where \f$\mathbf{u}_s\f$ is the surface value of the nonsliding SIA and 
\f$\mathbf{v}\f$ is the depth-independent velocity which goes in the 
SSA stress balance.  Note \f$\mathbf{u}_s\f$ is directly computable from
the temperature and geometry.  Recall \f$ f(|\mathbf{v}|) = 1 - (2/\pi)
\arctan\left(|\mathbf{v}|^2 v_0^{-2}\right) \f$ where \f$v_0=100\f$ m/a.  The 
unknown in the equation above is \f$\mathbf{v}\f$.

Our approach is to find \f$x = |\mathbf{v}|^2\f$ by solving a transcendental
equation by numerical root-finding.  Define 
\f$F(x) = (2/\pi) \arctan\left(x\,v_0^{-2}\right)\f$
so \f$f(|\mathbf{v}|) = 1 - F(x)\f$.  The scalar \f$x\f$ solves the scalar equation
  \f[ |\mathbf{U}_s - (1-F(x))\mathbf{u}_s|^2 = F(x)^2 x. \f]
We treat this as a root-finding problem \f$G(x)=0\f$ for \f$x\f$, where, with 
a little more rewriting,
  \f[ G(x) = (x - |\mathbf{u}_s|^2) F(x)^2 
             - 2 F(x) (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{u}_s
             - |\mathbf{U}_s - \mathbf{u}_s|^2. \f]
Note \f$F(0)=0\f$ and \f$F(x)\to 1\f$ as \f$x\to\infty\f$.  Thus \f$G(0)\le 0\f$ 
and \f$G(x) \sim x\f$ as \f$x\to\infty\f$.  Because \f$G(x)\f$ is continuous, we
know it has a nonnegative root.

We use the derivative
  \f[ G'(x) = F(x)^2 
              + 2 F'(x) \left[(x - |\mathbf{u}_s|^2) F(x) 
                              - (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{u}_s\right]\f]
in Newton's method
  \f[x_{n+1} = x_n - G(x_n)/G'(x_n).\f]

The method getGforInverse() actually computes \f$G(x)\f$ and 
\f$G'(x)\f$.  The initial guess is \f$x_0 = |\mathbf{u}_s|^2 
+ 4 (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{U}_s \f$ because that value
solves the equation that results from replacing the \f$F(x)\f$ term in the equation
\f$G(x)=0\f$ above by 1/2.  That is, \f$x_0\f$ would be the root of \f$G(x)\f$ if
\f$F(x)\f$ were equal to 1/2.

We iterate Newton's method to \f$10^{-12}\f$ relative tolerance and warn on nonconvergence.
 */
PetscErrorCode IceModel::computeFofVforInverse() {
  PetscErrorCode ierr;
  PetscScalar **us, **vs, **usSIA, **vsSIA, **imask, **H, **f;

  ierr = inv.usIn->get_array(us);    CHKERRQ(ierr);
  ierr = inv.vsIn->get_array(vs);    CHKERRQ(ierr);
  ierr = inv.usSIA->get_array(usSIA); CHKERRQ(ierr);
  ierr = inv.vsSIA->get_array(vsSIA); CHKERRQ(ierr);
  ierr = inv.invMask->get_array(imask); CHKERRQ(ierr);
  ierr = vH.get_array(H);    CHKERRQ(ierr);
  ierr = inv.fofv->get_array(f);     CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar USIAsqr = PetscSqr(usSIA[i][j]) + PetscSqr(vsSIA[i][j]);
      if (H[i][j] < 10.0) {
        f[i][j] = 1.0;  // treat all observed velocity as SIA if essentially no ice
      } else if (USIAsqr < PetscSqr(1.0/secpera)) {
        f[i][j] = 0.0;  // treat all observed velocity as sliding if ice is present
                        //   but there is almost no surface velocity
      } else if (imask[i][j] < 0.5) {
        f[i][j] = 1.0;  // treat all as SIA, but should not matter since SSA not applied there
      } else {
        // Newton's method on norm equation from convex combination using f
        const PetscScalar
          Usmus_x    = us[i][j] - usSIA[i][j],
          Usmus_y    = vs[i][j] - vsSIA[i][j],
          Usmus2     = PetscSqr(Usmus_x)     + PetscSqr(Usmus_y),
          Usmusdotus = Usmus_x * usSIA[i][j] + Usmus_y * vsSIA[i][j],
          UsmusdotUs = Usmus_x * us[i][j]    + Usmus_y * vs[i][j];
        PetscScalar x = USIAsqr + 4.0 * UsmusdotUs;  // initial guess solves nearby eqn
        PetscScalar G, Gprime;
        for (PetscInt k=0; k<4; ++k) {  // a fixed # of steps
          ierr = getGforInverse(x, Usmus2, Usmusdotus, USIAsqr, G, Gprime); CHKERRQ(ierr);
          x = x - G / Gprime;  // Newton; practicing unprotected steps here ...
        }

        // at this point x contains  x=|v|^2
        const PetscScalar  v0sqr = PetscSqr(100.0 / secpera);
        PetscScalar        fofv  = 1.0 - (2.0 / pi) * atan(x / v0sqr);

        // ad hoc: make SIA matter slightly less:
        //fofv = 0.9 * fofv;

        if (fofv < 0.0)    fofv = 0.0;
        if (fofv > 1.0)    fofv = 1.0;
        f[i][j] = fofv;
      }
    }
  }
  ierr = inv.fofv->end_access(); CHKERRQ(ierr);
  ierr = inv.usIn->end_access(); CHKERRQ(ierr);
  ierr = inv.vsIn->end_access(); CHKERRQ(ierr);
  ierr = inv.usSIA->end_access(); CHKERRQ(ierr);
  ierr = inv.vsSIA->end_access(); CHKERRQ(ierr);
  ierr = vH.end_access();    CHKERRQ(ierr);
  ierr = inv.invMask->end_access(); CHKERRQ(ierr);

  return 0;
}


//! Remove the SIA surface velocity from ice surface velocity to give depth-independent, membrane stress-balanced (SSA) velocity.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

This procedure computes a z-independent horizontal ice velocity, conceptually
a basal sliding velocity, from a specified surface horizontal ice velocity
(map).  The portion already found from SIA equations is removed.
 */
PetscErrorCode IceModel::removeSIApart() {
  PetscErrorCode ierr;

  // max amount before div by (1-f) will be bypassed:
  const PetscScalar maxfofv = 0.999;

  PetscScalar **us, **vs, **ub, **vb, **usSIA, **vsSIA, **fofv;
  ierr = inv.usIn->get_array(us);    CHKERRQ(ierr);
  ierr = inv.vsIn->get_array(vs);    CHKERRQ(ierr);
  ierr = inv.usSIA->get_array(usSIA); CHKERRQ(ierr);
  ierr = inv.vsSIA->get_array(vsSIA); CHKERRQ(ierr);
  ierr = inv.fofv->get_array(fofv);  CHKERRQ(ierr);
  ierr = vubarSSA.get_array(ub);    CHKERRQ(ierr);
  ierr = vvbarSSA.get_array(vb);    CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar f = fofv[i][j];
      if (f < maxfofv) {
        ub[i][j] = (us[i][j] - f * usSIA[i][j]) / (1.0 - f);
        vb[i][j] = (vs[i][j] - f * vsSIA[i][j]) / (1.0 - f);
      } else {
        // if f is 1, just treat all observed surface vel as SIA, so no sliding
        ub[i][j] = 0.0;
        vb[i][j] = 0.0;
      }
    }
  }
  ierr = inv.usIn->end_access(); CHKERRQ(ierr);
  ierr = inv.vsIn->end_access(); CHKERRQ(ierr);
  ierr = inv.usSIA->end_access(); CHKERRQ(ierr);
  ierr = inv.vsSIA->end_access(); CHKERRQ(ierr);
  ierr = inv.fofv->end_access(); CHKERRQ(ierr);
  ierr = vubarSSA.end_access(); CHKERRQ(ierr);
  ierr = vvbarSSA.end_access(); CHKERRQ(ierr);

  // SSA differential op will be applied; communicate
  ierr = vubarSSA.beginGhostComm(); CHKERRQ(ierr);
  ierr = vvbarSSA.beginGhostComm(); CHKERRQ(ierr);
  ierr = vubarSSA.endGhostComm(); CHKERRQ(ierr);
  ierr = vvbarSSA.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


//! Compute N, the effective pressure on the till, following basal model conventions.
PetscErrorCode IceModel::getEffectivePressureForInverse() {
  PetscErrorCode ierr;

  PetscScalar **N, **H, **Hmelt;
  ierr = inv.effPressureN->get_array(N);  CHKERRQ(ierr);
  ierr = vH.get_array(H);  CHKERRQ(ierr);
  ierr = vHmelt.get_array(Hmelt);  CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      N[i][j] = getEffectivePressureOnTill(H[i][j], Hmelt[i][j]); // iMbasal.cc
    }
  }
  ierr = vH.end_access();  CHKERRQ(ierr);
  ierr = vHmelt.end_access();  CHKERRQ(ierr);
  ierr = inv.effPressureN->end_access();  CHKERRQ(ierr);
  return 0;
}


//! Compute normalized velocity appropriate to inverse model.
/*!
Computes
  \f[ V = \frac{U}{|U|^{1-q} U_{th}^q} \f]
where \f$q=\f$ \c basal->pseudo_q and \f$U_{th}=\f$ \c basal->pseudo_u_threshold.
But regularizes the denominator by the replacement
  \f[ |U| \to \left(|U|^2 + (0.01 \text{ m/a})^2\right)^{1/2}. \f]
   
Used by computeTFAFromBasalShear() and computeTFAFromBasalShearNoReg(); see
comments for those methods.
 */
PetscErrorCode IceModel::getVfromUforInverse(
                  const PetscScalar U_x, const PetscScalar U_y,
                  PetscScalar &V_x, PetscScalar &V_y, PetscScalar &magVsqr) {
  const PetscScalar
              q       = basal->pseudo_q,
              Uth     = basal->pseudo_u_threshold,
              delta2  = PetscSqr(0.01/secpera),
              regUsqr = delta2 + PetscSqr(U_x)+ PetscSqr(U_y),
              denom   = pow(regUsqr, (1.0 - q) / 2.0) * pow(Uth,q);
  V_x = U_x / denom;
  V_y = U_y / denom;
  magVsqr = PetscSqr(V_x) + PetscSqr(V_y);
  return 0;
}


//! Compute the till friction angle for the given basal shear stress by minimizing an un-regularized objective.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

Normally computeTFAFromBasalShear() should be called to do this job with a
regularization scheme; see comments for that method.

Here we do something \e very simple, namely to divide magnitudes 
to get the till friction angle:
	\f[ \phi = \arctan\left(\frac{|\tau_b|}{N |\mathbf{V}|}\right).  \f]
This formula applies in \f$\Omega_O\f$; elsewhere we set \f$\phi=\phi_i\f$.
This calculation is only done where there is stencil width so that the SSA differential
operator could actually comput the basal shear stress \f$\tau_b\f$.

Finally we enforce the limits \c phi_min, \c phi_max.
 */
PetscErrorCode IceModel::computeTFAFromBasalShearNoReg(
                           const PetscScalar phi_low, const PetscScalar phi_high) {               
  PetscErrorCode ierr;
  PetscScalar    magVsqr, junk1, junk2;

  PetscScalar **phi, **oldphi, **ub, **vb, **imask, **N, **fofv, **taubx, **tauby;
  ierr = vtillphi.get_array(phi);  CHKERRQ(ierr);
  ierr = vubarSSA.get_array(ub);  CHKERRQ(ierr);
  ierr = vvbarSSA.get_array(vb);  CHKERRQ(ierr);
  ierr = inv.oldtillphi->get_array(oldphi);  CHKERRQ(ierr);
  ierr = inv.invMask->get_array(imask);  CHKERRQ(ierr);
  ierr = inv.effPressureN->get_array(N);  CHKERRQ(ierr);
  ierr = inv.fofv->get_array(fofv);  CHKERRQ(ierr);
  ierr = inv.taubxComputed->get_array(taubx); CHKERRQ(ierr);
  ierr = inv.taubyComputed->get_array(tauby); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (imask[i][j] < 1.5) {
        phi[i][j] = oldphi[i][j];  // no observed velocities or no stencil width
      } else {
        ierr = getVfromUforInverse(ub[i][j], vb[i][j], junk1, junk2, magVsqr); CHKERRQ(ierr);
        const PetscScalar
            magV    = sqrt(magVsqr),
            magtaub = sqrt(PetscSqr(taubx[i][j]) + PetscSqr(tauby[i][j]));
        // note   tau_b = - mu N V   while   mu = tan(phi)
        const PetscScalar idealphi = (180.0/pi) * atan(magtaub / (N[i][j] * magV));
        if (fofv[i][j] > 0.8) {  // if mostly SIA, gently turn off result and use old
          const PetscScalar lambda = (fofv[i][j] - 0.8) / 0.2;  // in [0,1]
          phi[i][j] =   (1.0 - lambda) * idealphi + lambda * oldphi[i][j];
        } else {
          phi[i][j] = idealphi;
        }
      }
      if (phi[i][j] > phi_high)  phi[i][j] = phi_high;
      if (phi[i][j] < phi_low)   phi[i][j] = phi_low;
    }
  }
  ierr = vtillphi.end_access();  CHKERRQ(ierr);
  ierr = vubarSSA.end_access();  CHKERRQ(ierr);
  ierr = vvbarSSA.end_access();  CHKERRQ(ierr);
  ierr = inv.oldtillphi->end_access();  CHKERRQ(ierr);
  ierr = inv.fofv->end_access();  CHKERRQ(ierr);
  ierr = inv.invMask->end_access();  CHKERRQ(ierr);
  ierr = inv.effPressureN->end_access();  CHKERRQ(ierr);
  ierr = inv.taubxComputed->end_access(); CHKERRQ(ierr);
  ierr = inv.taubyComputed->end_access(); CHKERRQ(ierr);

  return 0;
}

