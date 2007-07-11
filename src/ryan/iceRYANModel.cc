// Copyright (C) 2007 Ryan Woodard

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
#include <acml.h>
#include <acml_mv.h>
#include "grid.hh"
#include "materials.hh"
#include "iceEISModel.hh"

#include "iceRYANModel.hh"

IceRYANModel::IceRYANModel(IceGrid &g, IceType &i) 
  : IceEISModel(g,i) {  // do nothing; note derived classes must have constructors
}


PetscErrorCode IceRYANModel::setFromOptions() {
  PetscErrorCode ierr;
  char           temp, eisIIexpername[20], accname[20];
  PetscTruth     eisIIchosen, accchosen;

  ierr = IceEISModel::setFromOptions(); CHKERRQ(ierr);
  
  ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", eisIIexpername, 1, &eisIIchosen); CHKERRQ(ierr);
  temp = eisIIexpername[0];
  if ((eisIIchosen == PETSC_FALSE) || ((temp != 'h') && (temp != 'H'))) {
    SETERRQ(3,"pryan must be run with option '-eisII H'\n");
  }
  
  ierr = PetscOptionsGetString(PETSC_NULL, "-acc", accname, 1, &accchosen); CHKERRQ(ierr);
  if (accchosen == PETSC_TRUE) {
      temp = accname[0];
      if ((temp >= 'a') && (temp <= 'z'))   temp += 'A'-'a';  // capitalize if lower
      if ((temp >= 'A') && (temp <= 'F')) {
	     setAccName(temp);
      } else {
	     SETERRQ(2,"option -acc must have value A, B, C, D, E or F.\n");
      }
  } else {
    setAccName('A');  // default
    //SETERRQ(1,"option -acc must have a value\n");
  }
  return 0;
}


PetscErrorCode IceRYANModel::initFromOptions() {
  PetscErrorCode ierr;

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);

  // Ryan put this at the end of  IceEISModel::initAccumTs()
  setM_max_g(M_max);
  setS_b_g(S_b);
  setR_el_g(R_el);
  
  ierr = initRandomnessACML(); CHKERRQ(ierr);

  return 0;
}



PetscErrorCode IceRYANModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscInt tempskipCount, const char adaptReason,
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {
  
  PetscErrorCode ierr;
  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,"       YEAR    VOL           AREA\n");  CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, "%11.3f %13.8f %12.7f",
		      year, volume_kmcube/1.0e6, area_kmsquare/1.0e6); CHKERRQ(ierr);
  }
  return 0;
}


void IceRYANModel::setExperName_g(char name) {
  expername_g = name;
}

char IceRYANModel::getExperName_g() {
  return expername_g;
}

void IceRYANModel::setAccName(char name) {
  accname = name;
}

char IceRYANModel::getAccName() {
  return accname;
}

void IceRYANModel::setM_max_g(PetscScalar M_max) {
  M_max_g = M_max;
}

void IceRYANModel::setS_b_g(PetscScalar S_b) {
  S_b_g = S_b;
}

void IceRYANModel::setR_el_g(PetscScalar R_el) {
  R_el_g = R_el;
}

PetscScalar IceRYANModel::getM_max_g() {
  return M_max_g;
}

PetscScalar IceRYANModel::getS_b_g() {
  return S_b_g;
}

PetscScalar IceRYANModel::getR_el_g() {
  return R_el_g;
}


PetscErrorCode IceRYANModel::initRandomnessACML() {
  PetscErrorCode ierr;
  PetscInt       *seed, *state;
  PetscScalar    gauss_xmu, gauss_var, lognormal_xmu, lognormal_var;

  randomnessacml.lseed = 1;
  randomnessacml.lstate = 641;

  // original M_max = 0.5 or 406 kg m^2 a^-1
  // mean of positive values of Arthern/Vaughan's
  // accum. data is 185 kg m^2 a^-1
  gauss_xmu = 0.5 * 185. / 406.;
  gauss_var = pow(gauss_xmu / 3.0, 2);
  
  // xmu and var from normal fit to 
  // log of positive values of Arthern/Vaughan's data
  //lognormal_xmu = 4.858961;
  lognormal_xmu = -1.8405;
  lognormal_var = pow(1.227035, 2);

  switch (getAccName()) {
    case 'B':
      randomnessacml.n = grid.p->Mx * grid.p->My; // 1782225;
      randomnessacml.xmu = gauss_xmu;
      randomnessacml.var = gauss_var;
      break;
    case 'C':
      randomnessacml.n = 1; // recomputed at each time step
                            // since we do not know how many
                            // time steps there will be (adaptive dt)
      randomnessacml.xmu = gauss_xmu;
      randomnessacml.var = gauss_var;
      break;
    case 'D':
      randomnessacml.n = grid.p->Mx * grid.p->My; // 1782225;
      randomnessacml.xmu = lognormal_xmu;
      randomnessacml.var = lognormal_var;
      break;
    case 'E':
      randomnessacml.n = 1; // recomputed at each time step
                            // since we do not know how many
                            // time steps there will be (adaptive dt)
      randomnessacml.xmu = lognormal_xmu;
      randomnessacml.var = lognormal_var;
      break;
    default: // This line should never be reached.
      break;
  }
  
  PetscMalloc(randomnessacml.lseed * sizeof(PetscInt), &seed);
  ISCreateGeneral(PETSC_COMM_SELF, randomnessacml.lseed, seed, &randomnessacml.seed);
  // Note that ISCreateGeneral() has made a copy of seed
  // so we may (and generally should) free seed
  PetscFree(seed);

  PetscMalloc(randomnessacml.lstate * sizeof(PetscInt), &state);
  ISCreateGeneral(PETSC_COMM_SELF, randomnessacml.lstate, state, &randomnessacml.state);
  PetscFree(state);

  ierr = VecCreateSeq(PETSC_COMM_SELF, randomnessacml.n, &randomnessacml.x);

  randomnessacml.genid = 3;
  randomnessacml.subid = 1;

  ierr = ISGetIndices(randomnessacml.seed, &seed); CHKERRQ(ierr);
  ierr = ISGetIndices(randomnessacml.state, &state); CHKERRQ(ierr);

  PetscScalar mySeed; // Always want to spread it.
  PetscTruth  paramSet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ACMLseed", &mySeed, &paramSet); CHKERRQ(ierr);
  if (paramSet == PETSC_TRUE)
    seed[0] = (PetscInt) mySeed;
  else
    seed[0] = 17;
  
  drandinitialize(randomnessacml.genid, randomnessacml.subid,
		  seed, &randomnessacml.lseed,
		  state, &randomnessacml.lstate,
		  &randomnessacml.info);

  ierr = ISRestoreIndices(randomnessacml.seed, &seed); CHKERRQ(ierr);
  ierr = ISRestoreIndices(randomnessacml.state, &state); CHKERRQ(ierr);

  if (randomnessacml.info != 0) {
    SETERRQ(1,"Error in drandinitialize:  info != 0\n");
  }
  return 0;
}


PetscErrorCode IceRYANModel::additionalAtStartTimestep() {
  PetscErrorCode    ierr;
  PetscScalar       M_max, S_b;
  PetscScalar       **accum, *x;
  PetscInt          *state, k;

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = ISGetIndices(randomnessacml.state, &state); CHKERRQ(ierr);
  ierr = VecGetArray(randomnessacml.x, &x); CHKERRQ(ierr);

  const PetscScalar R_el = getR_el_g();

  switch (getAccName()) {
    case 'B':
    case 'C':
      //printf("Normal Mmax acc.\n");
      drandgaussian(randomnessacml.n, 
		    randomnessacml.xmu, randomnessacml.var,
		    state, x, &randomnessacml.info);
      break;
    case 'D':
    case 'E':
      //printf("Lognormal Mmax acc.\n");
      drandlognormal(randomnessacml.n, 
		     randomnessacml.xmu, randomnessacml.var,
		     state, x, &randomnessacml.info);
      break;
  }
  if (randomnessacml.info != 0) {
    SETERRQ(1,"Error in drandinitialize:  info != 0\n");
  }
  switch (getAccName()) {
    case 'B':
    case 'D':
      k = 0;
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
	for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	  // r is distance from center of grid
	  const PetscScalar r = sqrt( PetscSqr(-grid.p->Lx + grid.p->dx*i)
				      + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
	  if (r < R_el) {
	    // set accumulation
	    //accum[i][j] = 1.58444e-08;
	    accum[i][j] = x[k] / secpera;
	    k++;
// 	    printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",
// 		   getM_max_g(), x[k], M_max, S_b, accum[i][j]);
	  }
	}
      }
      break;
    case 'C':
    case 'E':
      M_max = x[0] / secpera;
      if (M_max < 0.0) M_max = 0.0;
      S_b = getS_b_g() * M_max / getM_max_g();
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
	for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	  // r is distance from center of grid
	  const PetscScalar r = sqrt( PetscSqr(-grid.p->Lx + grid.p->dx*i)
				      + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
	  // set accumulation
	  //accum[i][j] = 1.58444e-08;
	  accum[i][j] = PetscMin(M_max, S_b * (R_el - r));  // formula (7) in (Payne et al 2000)
	  //printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n", getM_max_g(), x[0], M_max, S_b, accum[i][j]);
	}
      }
      break;
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = ISRestoreIndices(randomnessacml.state, &state); CHKERRQ(ierr);
  ierr = VecRestoreArray(randomnessacml.x, &x); CHKERRQ(ierr);
  return 0;
}

