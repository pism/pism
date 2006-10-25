// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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

// Process command line options.  This should be called by a driver
// program, if it would like to use command line options, _after_
// all driver overrides are performed.

PetscErrorCode
IceModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscScalar my_maxdt, my_mu, my_startYear, my_runYears, my_endYear;
  PetscScalar macRTol, my_nu, maceps;
  PetscTruth my_useMacayealVelocity, my_useConstantNu, macRTolSet, macepsSet,
             maxdtSet, startYearSet, runYearsSet, endYearSet, verbose,
             noMassBal, noTemp, bedDefiso, bedDef, bedDeflc, isoflux, muSet, 
             nospokesSet, tempskipSet, oceanKillSet;
  PetscInt nospokeslevel, ts;

  ierr = PetscOptionsGetScalar(PETSC_NULL, "-adapt_ratio", &adaptTimeStepRatio,
                               PETSC_NULL); CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL, "-bed_def_lc", &bedDeflc); CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL, "-bed_def", &bedDef); CHKERRQ(ierr);
  if (bedDef == PETSC_TRUE) {   bedDeflc = PETSC_TRUE; }
  
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
  
  ierr = PetscOptionsHasName(PETSC_NULL, "-no_mass_bal", &noMassBal); CHKERRQ(ierr);
  if (noMassBal == PETSC_TRUE) {
    setDoMassBal(PETSC_FALSE);
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

// note "-of" is in use for output file format

// note "-regrid" is in use for regrid file name

// note "-regrid_vars" is in use for regrid variable names

  /* This controls how many mass balance steps per temp step */
  ierr = PetscOptionsGetInt(PETSC_NULL, "-tempskip", &ts, &tempskipSet); CHKERRQ(ierr);
  if (tempskipSet == PETSC_TRUE)
    setTempskip(ts);
  
  // verbose: the summary at each time step is more complete (see summary() in 
  //     iMutil.cc)
  ierr = PetscOptionsHasName(PETSC_NULL, "-verbose", &verbose); CHKERRQ(ierr);
  if (verbose == PETSC_TRUE) {
    setBeVerbose(PETSC_TRUE);
  }

  // Run length options
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ys", &my_startYear, &startYearSet); CHKERRQ(ierr);
  if (startYearSet == PETSC_TRUE) {
    setStartYear(my_startYear);
  }
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-y", &my_runYears, &runYearsSet); CHKERRQ(ierr);
  if (runYearsSet == PETSC_TRUE) {
    ierr = setRunYears(my_runYears); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ye", &my_endYear, &endYearSet); CHKERRQ(ierr);
  if (endYearSet == PETSC_TRUE) {
    ierr = setEndYear(my_endYear); CHKERRQ(ierr);
  }
    
  return 0;
}
