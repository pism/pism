// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include "iceModel.hh"

// Process command line options.  This is called by a driver
// program, assuming it would like to use command line options.

/* Note no options to directly set dx, dy, dz, year as the user should not directly set these quantities; 
see options for Mx,My,Mz,Mbz,Lx,Ly,Lz,ys,ye,y. */

PetscErrorCode  IceModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscScalar my_maxdt, my_mu;
  PetscScalar macRTol, my_nu, maceps, regVelSchoof, regLengthSchoof;
  PetscTruth my_useMacayealVelocity, my_useConstantNu, macRTolSet, macepsSet,
             maxdtSet, superpose,
             noMassConserve, noTemp, bedDefiso, bedDeflc, isoflux, muSet, 
             nospokesSet, oceanKillSet, tempskipSet, regVelSchoofSet, regLengthSchoofSet,
             MxSet, MySet, MzSet, MbzSet;
  PetscInt nospokeslevel, my_Mx, my_My, my_Mz, my_Mbz;

  // OptionsBegin/End probably has no effect for now, but perhaps some day PETSc will show a GUI which
  // allows users to set options using this.
  ierr = PetscOptionsBegin(grid.com,PETSC_NULL,"IceModel options (in PISM)",PETSC_NULL); CHKERRQ(ierr);

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-adapt_ratio", &adaptTimeStepRatio,
                               PETSC_NULL); CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL, "-bed_def_lc", &bedDeflc); CHKERRQ(ierr);  
  ierr = PetscOptionsHasName(PETSC_NULL, "-bed_def_iso", &bedDefiso); CHKERRQ(ierr);
  if (bedDefiso == PETSC_TRUE) {
    setDoBedIso(PETSC_TRUE);
  }
  if ((bedDefiso == PETSC_TRUE) || (bedDeflc == PETSC_TRUE)) {
    setDoBedDef(PETSC_TRUE);
  }

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-constant_nu", &my_nu, &my_useConstantNu); CHKERRQ(ierr);
  // user gives nu in MPa yr (e.g. Ritz value is 30.0)
  if (my_useConstantNu == PETSC_TRUE) {
    setConstantNuForMacAyeal(my_nu  * 1.0e6 * secpera);
  }

  // regular size viewers
  ierr = PetscOptionsGetString(PETSC_NULL, "-d", diagnostic, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  if (showViewers == PETSC_FALSE) {
    strcpy(diagnostic, "\0");
  }

  // big viewers (which will override regular viewers)
  ierr = PetscOptionsGetString(PETSC_NULL, "-dbig", diagnosticBIG, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  if (showViewers == PETSC_FALSE) {
    strcpy(diagnosticBIG, "\0");
  }

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-e", &enhancementFactor, PETSC_NULL); CHKERRQ(ierr);

// note "-gk" is in use for specifying Goldsby-Kohlstedt ice

// note "-id" is in use for sounding location

// note "-if" is in use for input file name

  // This switch turns off vertical integration in the isothermal case.  That
  // is, the horizontal flux of ice is computed as an analytical function of the
  // thickness and the surface slope.  The Glen power n=3 and a fixed softness
  // parameter A = 10^{-16} Pa^{-3} a^{-1} are used.  These are set in
  // IceModel::setDefaults().
  ierr = PetscOptionsHasName(PETSC_NULL, "-isoflux", &isoflux); CHKERRQ(ierr);
  setIsothermalFlux(isoflux);

// note "-jd" is in use for sounding location

// note "-kd" is in use for horizontal slicing (in viewers and dumpToFileMatlab)

// note -Lx, -Ly, -Lz are all checked in [iMutil.cc]IceModel::afterInitHook()

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-maxdt", &my_maxdt, &maxdtSet); CHKERRQ(ierr);
  if (maxdtSet == PETSC_TRUE) {
    setMaxTimeStepYears(my_maxdt);
  }

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-mu_sliding", &my_mu, &muSet); CHKERRQ(ierr);
  if (muSet == PETSC_TRUE) {
    setMuSliding(my_mu);
  }

  ierr = PetscOptionsHasName(PETSC_NULL, "-mv", &my_useMacayealVelocity); CHKERRQ(ierr);
  if (my_useMacayealVelocity == PETSC_TRUE) {
    setUseMacayealVelocity(my_useMacayealVelocity);
  }
  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-mv_eps", &maceps, &macepsSet); CHKERRQ(ierr);
  if (macepsSet == PETSC_TRUE) {
    setMacayealEpsilon(maceps);
  }
  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-mv_rtol", &macRTol, &macRTolSet); CHKERRQ(ierr);
  if (macRTolSet == PETSC_TRUE) {
    setMacayealRelativeTolerance(macRTol);
  }
  
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mx", &my_Mx, &MxSet); CHKERRQ(ierr);
  if (MxSet == PETSC_TRUE)   grid.p->Mx = my_Mx;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-My", &my_My, &MySet); CHKERRQ(ierr);
  if (MySet == PETSC_TRUE)   grid.p->My = my_My;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mz", &my_Mz, &MzSet); CHKERRQ(ierr);
  if (MzSet == PETSC_TRUE)   grid.p->Mz = my_Mz;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mbz", &my_Mbz, &MbzSet); CHKERRQ(ierr);
  if (MbzSet == PETSC_TRUE)   grid.p->Mbz = my_Mbz;

  ierr = PetscOptionsHasName(PETSC_NULL, "-no_mass", &noMassConserve); CHKERRQ(ierr);
  if (noMassConserve == PETSC_TRUE) {
    setDoMassConserve(PETSC_FALSE);
  }

  // -no_spokes K for K=0,1,2,... turns on smoothing of spokes by smoothing Sigma (e.g. in EISMINT experiment F)
  // values K>3 not recommended (lots of communication!)
  ierr = PetscOptionsGetInt(PETSC_NULL, "-no_spokes", &nospokeslevel, &nospokesSet); CHKERRQ(ierr);
  if (nospokesSet == PETSC_TRUE)
    setNoSpokes(nospokeslevel);

  ierr = PetscOptionsHasName(PETSC_NULL, "-no_temp", &noTemp); CHKERRQ(ierr);
  if (noTemp == PETSC_TRUE) {
    setDoTemp(PETSC_FALSE);
  }

// note "-o" is in use for output file name

  // whether or not to kill ice if original condition was ice-free ocean
  ierr = PetscOptionsHasName(PETSC_NULL, "-ocean_kill", &oceanKillSet); CHKERRQ(ierr);
  if (oceanKillSet == PETSC_TRUE) {
    setOceanKill(PETSC_TRUE);
  }

// note "-of" is in use for output file format; see iMIO.cc

  // use a plastic basal till mechanical model
  ierr = PetscOptionsHasName(PETSC_NULL, "-plastic", &doPlasticTill); CHKERRQ(ierr);

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-reg_vel_schoof", &regVelSchoof, &regVelSchoofSet); CHKERRQ(ierr);
  if (regVelSchoofSet == PETSC_TRUE) {
    setRegularizingVelocitySchoof(regVelSchoof/secpera);
  }
  
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-reg_length_schoof", &regLengthSchoof, &regLengthSchoofSet); CHKERRQ(ierr);
  if (regLengthSchoofSet == PETSC_TRUE) {
    setRegularizingLengthSchoof(regLengthSchoof);
  }
  
// note "-regrid" is in use for regrid file name; see iMregrid.cc

// note "-regrid_vars" is in use for regrid variable names; see iMregrid.cc

  // apply "glaciological superposition to low order", i.e. add SIA results to those of MacAyeal equations
  // where DRAGGING
  ierr = PetscOptionsHasName(PETSC_NULL, "-super", &superpose); CHKERRQ(ierr);
  setDoSuperpose(superpose);

  /* This controls allows more than one mass continuity steps per temperature/age step */
  ierr = PetscOptionsGetInt(PETSC_NULL, "-tempskip", &tempskipMax, &tempskipSet); CHKERRQ(ierr);
  if (tempskipSet == PETSC_TRUE) {
    doTempSkip = PETSC_TRUE;
  }

  // verbosity options: more info to standard out.  see iMutil.cc
  ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

  // -ys, -ye, -y options read in setStartRunEndYearsFromOptions()
 
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  return 0;
}
