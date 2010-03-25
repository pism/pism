// Copyright (C) 2006-2010 Ed Bueler, Constantine Khroulev, and Jed Brown
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

#include <petsc.h>
#include "pism_const.hh"
#include "grid.hh"
#include "iceModelVec.hh"

#include "stressBalanceSSA.hh"


newSSAStrengthExtension::newSSAStrengthExtension() {
  min_thickness = 50.0;   // m
          // minimum thickness (for SSA velocity computation) at which 
          // NuH switches from vertical integral to constant value
          // this value strongly related to calving front
          // force balance, but the geometry itself is not affected by this value
  const PetscReal
    DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8,  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
    DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / secpera) / (100.0 * 1.0e3);  // typical strain rate is 100 m/yr per 
  nuH = min_thickness * DEFAULT_CONSTANT_HARDNESS_FOR_SSA
                       / (2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); // Pa s m
          // COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
          //          30 MPa yr for \bar\nu
}


PetscErrorCode newSSAStrengthExtension::set_notional_strength(PetscReal my_nuH) {
  nuH = my_nuH;
  return 0;
}

PetscErrorCode newSSAStrengthExtension::set_min_thickness(PetscReal my_min_thickness) {
  min_thickness = my_min_thickness;
  return 0;
}


StressBalanceSSA::StressBalanceSSA(
           IceGrid* g, IceFlowLaw* ssa_ice, IceBasalResistancePlasticLaw* ssa_basal,
           IceModelVec2 *ssa_tauc, IceModelVec2 *ssa_mask, IceModelVec2 *ssa_hardav) {
  initAndAllocate(g);
  ice = ssa_ice;
  basal = ssa_basal;
  vtauc = ssa_tauc;
  vmask = ssa_mask;
  vhardav = ssa_hardav;
}


StressBalanceSSA::~StressBalanceSSA() {
  deallocate();
}


PetscErrorCode StressBalanceSSA::initAndAllocate(IceGrid* g) {
  PetscErrorCode ierr;
  
  grid = g;

  config.init("pism_config", grid->com, grid->rank);
  char alt_config[PETSC_MAX_PATH_LEN];
  PetscTruth use_alt_config;
  ierr = PetscOptionsGetString(PETSC_NULL, "-config", alt_config, 
                               PETSC_MAX_PATH_LEN, &use_alt_config); CHKERRQ(ierr);
  if (use_alt_config) {
    ierr = config.read(alt_config); CHKERRQ(ierr);
  } else {
    ierr = config.read(PISM_DefaultConfigFile); CHKERRQ(ierr);
  }

  char overridename[PETSC_MAX_PATH_LEN];
  PetscTruth use_override_config;
  ierr = PetscOptionsGetString(PETSC_NULL, "-config_override", overridename,
			       PETSC_MAX_PATH_LEN, &use_override_config);
  if (use_override_config) {
    NCConfigVariable overrides;
    overrides.init("pism_overrides", grid->com, grid->rank);
    ierr = overrides.read(overridename); CHKERRQ(ierr);
    config.import_from(overrides);
  }

  config.print();

  // setup (classical) SSA tools
  const PetscInt M = 2 * grid->Mx * grid->My;

  ierr = MatCreateMPIAIJ(grid->com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &SSAStiffnessMatrix); CHKERRQ(ierr);

  ierr = VecCreateMPI(grid->com, PETSC_DECIDE, M, &SSAX); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);

  ierr = VecCreateSeq(PETSC_COMM_SELF, M, &SSAXLocal);
  ierr = VecScatterCreate(SSAX, PETSC_NULL, SSAXLocal, PETSC_NULL,
                          &SSAScatterGlobalToLocal); CHKERRQ(ierr);

  ierr = KSPCreate(grid->com, &SSAKSP); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(SSAKSP); CHKERRQ(ierr);

  // initial guesses of SSA velocities
  ierr = ssavel.create(grid, "bar_ssa", true); // components are ubar_ssa and vbar_ssa
  ierr = ssavel.set_attrs("internal_restart", "SSA model ice velocity in the X direction",
                            "m s-1", "", 0); CHKERRQ(ierr);
  ierr = ssavel.set_attrs("internal_restart", "SSA model ice velocity in the Y direction",
                            "m s-1", "", 1); CHKERRQ(ierr);
  ierr = ssavel.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ierr = vaveragedhardness.create(*grid, "hardavSSA", true); CHKERRQ(ierr);
  const PetscScalar power = 1.0 / ice->exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);
  ierr = vaveragedhardness.set_attrs("ssa_internal",
                                     "vertical average of ice hardness",
			             unitstr, ""); CHKERRQ(ierr);
  ierr = vaveragedhardness.set_attr("valid_min",0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode StressBalanceSSA::deallocate() {
  PetscErrorCode ierr;

  ierr = KSPDestroy(SSAKSP); CHKERRQ(ierr);
  ierr = MatDestroy(SSAStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecDestroy(SSAX); CHKERRQ(ierr);
  ierr = VecDestroy(SSARHS); CHKERRQ(ierr);
  ierr = VecDestroy(SSAXLocal); CHKERRQ(ierr);
  ierr = VecScatterDestroy(SSAScatterGlobalToLocal); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode StressBalanceSSA::setGuessZero() {
  PetscErrorCode ierr;
  ierr = ssavel.set(0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode StressBalanceSSA::setGuess(IceModelVec2S *ubar_guess, IceModelVec2S *vbar_guess) {
  PetscErrorCode ierr;
  ierr = ssavel.set_component(0, *ubar_guess); CHKERRQ(ierr);
  ierr = ssavel.set_component(1, *vbar_guess); CHKERRQ(ierr);
  return 0;
}


PetscScalar StressBalanceSSA::basalDragx(PetscScalar **tauc, PISMVector2 **uv,
                                         PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], uv[i][j].u, uv[i][j].v);
}


PetscScalar StressBalanceSSA::basalDragy(PetscScalar **tauc, PISMVector2 **uv,
                                         PetscInt i, PetscInt j) const {                         
  return basal->drag(tauc[i][j], uv[i][j].u, uv[i][j].v);
};


/*!
Compute compoents of the basal stress applied to the base of the ice:
  \f[ \tau_{b,x} = - C(\tau_c,u,v) u, \f]
  \f[ \tau_{b,y} = - C(\tau_c,u,v) v, \f]
 */
PetscErrorCode StressBalanceSSA::getBasalStress(IceModelVec2 *vbs_x, IceModelVec2 *vbs_y) {
  PetscErrorCode ierr;
  
  PISMVector2 **uv;
  PetscScalar **tauc, **mask, **bs_x, **bs_y;
  ierr = ssavel.get_array(uv); CHKERRQ(ierr);
  ierr = vtauc->get_array(tauc); CHKERRQ(ierr);
  ierr = vmask->get_array(mask); CHKERRQ(ierr);
  ierr = vbs_x->get_array(bs_x); CHKERRQ(ierr);
  ierr = vbs_y->get_array(bs_y); CHKERRQ(ierr);

  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      bs_x[i][j] = - basalDragx(tauc, uv, i, j) * u[i][j];
      bs_y[i][j] = - basalDragy(tauc, uv, i, j) * v[i][j];
    }
  }

  ierr = ssavel.end_access(); CHKERRQ(ierr);
  ierr = vtauc->end_access(); CHKERRQ(ierr);
  ierr = vmask->end_access(); CHKERRQ(ierr);
  ierr = vbs_x->end_access(); CHKERRQ(ierr);
  ierr = vbs_y->end_access(); CHKERRQ(ierr);
  
  return 0;
}

