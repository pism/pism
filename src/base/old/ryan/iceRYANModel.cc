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
// #include <acml_mv.h>  only seems available under 64 bit, which seems incompatible with my machine
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../eismint/iceEISModel.hh"

#include "iceRYANModel.hh"

IceRYANModel::IceRYANModel(IceGrid &g, IceType &i) 
  : IceEISModel(g,i) {  // do nothing; note derived classes must have constructors
}


PetscErrorCode IceRYANModel::setFromOptions() {
  PetscErrorCode ierr;
  char           temp, eisIIexpername[20], longaccname[20];
  PetscTruth     eisIIchosen, accchosen;

  ierr = IceEISModel::setFromOptions(); CHKERRQ(ierr);
  
  // user may optionally give "-eisII H" but this is not necessary
  expername = 'H';  // set IceEISModel::expername directly
  ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", eisIIexpername, 1, &eisIIchosen); CHKERRQ(ierr);
  temp = eisIIexpername[0];
  if ((eisIIchosen == PETSC_TRUE) && ((temp != 'h') && (temp != 'H'))) {
    SETERRQ(1,"pryan only runs EISMINT II experiment H; other experiments not allowed\n");
  }

  // option "-acc ?" for perturbation modes ?=A,B,C,D,E
  // note "-acc A" *or* no given option "-acc" means plain vanilla EISMINT II H
  accname = 'A';  
  ierr = PetscOptionsGetString(PETSC_NULL, "-acc", longaccname, 1, &accchosen); CHKERRQ(ierr);
  if (accchosen == PETSC_TRUE) {
      temp = longaccname[0];
      if ((temp >= 'a') && (temp <= 'z'))   temp += 'A'-'a';  // capitalize if lower
      if ((temp >= 'A') && (temp <= 'E')) {
	     accname = temp;
      } else {
	     SETERRQ(2,"option -acc must have value A, B, C, D, or E.\n");
      }
  }

  mySeed = 17;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-ACMLseed", &mySeed, PETSC_NULL); CHKERRQ(ierr);

  ierr = verbPrintf(4,grid.com,"accname = %c and mySeed = %d in IceRYANModel\n",accname,mySeed);
          CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceRYANModel::initFromOptions() {
  PetscErrorCode ierr;

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);

  // set the parameters for IceEISModel just to be sure they have EISMINT II experiment H values
  S_b = 1e-2 * 1e-3 / secpera, // Grad of accum rate change
  M_max = 0.5 / secpera,       // Max accumulation (at center)
  R_el = 450e3;                // Distance to equil line (accum=0)
  
  ierr = initRandomnessACML(); CHKERRQ(ierr);

  ierr = verbPrintf(4,grid.com,"IceRYANModel initialized\n"); CHKERRQ(ierr);
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
  lognormal_xmu = -1.8405;
  lognormal_var = pow(1.227035, 2);

  switch (accname) {
    case 'A': // set up randomness same for A and B; just don't *use* it for A
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
    default:
      SETERRQ(1,"accname invalid");
      break;
  }
  
  PetscMalloc(randomnessacml.lseed * sizeof(PetscInt), &seed);
  ISCreateGeneral(PETSC_COMM_SELF, randomnessacml.lseed, seed, &randomnessacml.seed);
  // Note that ISCreateGeneral() has made a copy of seed
  // so we may (and generally should) free seed
  ierr = PetscFree(seed); CHKERRQ(ierr);

  PetscMalloc(randomnessacml.lstate * sizeof(PetscInt), &state);
  ISCreateGeneral(PETSC_COMM_SELF, randomnessacml.lstate, state, &randomnessacml.state);
  ierr = PetscFree(state); CHKERRQ(ierr);

  ierr = VecCreateSeq(PETSC_COMM_SELF, randomnessacml.n, &randomnessacml.x);

  randomnessacml.genid = 3;
  randomnessacml.subid = 1;

  ierr = ISGetIndices(randomnessacml.seed, &seed); CHKERRQ(ierr);
  ierr = ISGetIndices(randomnessacml.state, &state); CHKERRQ(ierr);

  seed[0] = mySeed;
  
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
  
  if (accname == 'A')  // do nothing additional; just continue EISMINT II experiment H
    return 0;

  PetscScalar       **accum, *x;
  PetscInt          *state, k;

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = ISGetIndices(randomnessacml.state, &state); CHKERRQ(ierr);
  ierr = VecGetArray(randomnessacml.x, &x); CHKERRQ(ierr);

  switch (accname) {
    case 'B':
    case 'C':
      drandgaussian(randomnessacml.n, 
		    randomnessacml.xmu, randomnessacml.var,
		    state, x, &randomnessacml.info);
      break;
    case 'D':
    case 'E':
      drandlognormal(randomnessacml.n, 
		     randomnessacml.xmu, randomnessacml.var,
		     state, x, &randomnessacml.info);
      break;
    default:
      SETERRQ(1,"accname invalid");
      break;
  }
  if (randomnessacml.info != 0) {
    SETERRQ(2,"Error in drandinitialize:  info != 0\n");
  }

  switch (accname) {
    case 'B':
    case 'D':
      k = 0;
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
          // r is distance from center of grid
          const PetscScalar r = sqrt(  PetscSqr(-grid.p->Lx + grid.p->dx*i)
                                     + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
          if (r < R_el) {
            // set accumulation
            accum[i][j] = x[k] / secpera;
            k++;
          }
        }
      }
      break;
    case 'C':
    case 'E':
      PetscScalar my_M_max, my_S_b;
      my_M_max = x[0] / secpera;  // apply randomness uniformly
      if (my_M_max < 0.0)  my_M_max = 0.0;
      my_S_b = S_b * my_M_max / M_max;
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
          // r is distance from center of grid
          const PetscScalar r = sqrt(  PetscSqr(-grid.p->Lx + grid.p->dx*i)
                                     + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
          // set accumulation
          accum[i][j] = PetscMin(my_M_max, my_S_b * (R_el - r));  // formula (7) in (Payne et al 2000)
        }
      }
      break;
    default:
      SETERRQ(3,"accname invalid");
      break;
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = ISRestoreIndices(randomnessacml.state, &state); CHKERRQ(ierr);
  ierr = VecRestoreArray(randomnessacml.x, &x); CHKERRQ(ierr);

  ierr = verbPrintf(4,grid.com,"additionalAtStartTimestep() in IceRYANModel done; "
                               "accumulation set for accname = %c\n",accname); 
            CHKERRQ(ierr);
  return 0;
}

