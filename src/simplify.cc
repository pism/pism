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

static char help[] =
  "Ice sheet driver for thermocoupled EISMINT II and other simplified geometry\n"
  "ice sheet simulations.\n";
// note this source defines derived class IceEISModel *and* provides driver itself

#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

class IceEISModel : public IceModel {
public:
    IceEISModel(IceGrid &g, IceType &i);
    virtual PetscErrorCode initFromOptions();
    void     setExperName(char);
    char     getExperName();
    void     setflowlawNumber(PetscInt);
    PetscInt getflowlawNumber();
    
private:
    char       expername;
    PetscTruth inFileSet;
    PetscInt   flowlawNumber;
 
    virtual PetscScalar basal(const PetscScalar x, const PetscScalar y,
         const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
         const PetscScalar mu);
    PetscErrorCode setExperNameFromOptions();
    PetscErrorCode applyDefaultsForExperiment();
    PetscErrorCode initAccumTs();
    PetscErrorCode fillintemps();
};


IceEISModel::IceEISModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {
}


void IceEISModel::setExperName(char name) {
  expername = name;
}


char IceEISModel::getExperName() {
  return expername;
}


void IceEISModel::setflowlawNumber(PetscInt law) {
  flowlawNumber = law;
}


PetscInt IceEISModel::getflowlawNumber() {
  return flowlawNumber;
}

// reimplement IceModel::basal:
PetscScalar IceEISModel::basal(const PetscScalar x, const PetscScalar y,
      const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
      const PetscScalar mu) {
  // note this version ignors mu
  
  if (getExperName() == '0') { // ISMIP-HEINO
    //PetscErrorCode  ierr = PetscPrintf(grid.com, 
    //        "   [IceEISModel::basal called with:   x=%f, y=%f, H=%f, T=%f, alpha=%f]\n",
    //        x,y,H,T,alpha);  CHKERRQ(ierr);
    if (T + ice.beta_CC_grad * H > DEFAULT_MIN_TEMP_FOR_SLIDING) {
      // set coords according to ISMIP-HEINO convention
      const PetscScalar  xIH = x + grid.p->Lx,
                         yIH = y + grid.p->Ly;
      if (   ( (xIH >= 2300e3-1) && (xIH <= 3300e3+1)
               && (yIH >= 1500e3-1) && (yIH <= 2500e3+1) ) 
          || ( (xIH >= 3300e3-1) && (xIH <= 4000e3+1)
               && (yIH >= 1900e3-1) && (yIH <= 2100e3+1) ) ) {
        // "soft sediment"; C_S = 500 a^-1
        return (500 / secpera) * H;
      } else {
        // "hard rock"; C_R = 10^5 a^-1; alpha = |grad h|
        return (1e5 /secpera) * H * alpha * alpha;
      }
    } else
      return 0.0;  // not at pressure melting
  } else
    return 0.0;  // zero sliding for other tests
}


PetscErrorCode IceEISModel::setExperNameFromOptions() {
  PetscErrorCode      ierr;
  char                temp, eisIIexpername[20], ismipexpername[20];
  PetscTruth          EISIIchosen, ISMIPchosen;

  /* This option determines the single character name of EISMINT II experiment:
  "-eisII F", for example. */
  ierr = PetscOptionsGetString(PETSC_NULL, "-eisII", eisIIexpername, 1, &EISIIchosen); CHKERRQ(ierr);
  /* This option chooses ISMIP; "-ismip H" is ISMIP-HEINO */
  ierr = PetscOptionsGetString(PETSC_NULL, "-ismip", ismipexpername, 1, &ISMIPchosen); CHKERRQ(ierr);
  
  if (EISIIchosen == PETSC_TRUE) {
    temp = eisIIexpername[0];
    if ((temp >= 'a') && (temp <= 'z'))   temp += 'A'-'a';  // capitalize if lower
    setExperName(temp);
  } else if (ISMIPchosen == PETSC_TRUE) {
    temp = ismipexpername[0];
    if ((temp == 'H') || (temp == 'h')) // ISMIP-HEINO
        setExperName('0');
    else if ((temp == 'O') || (temp == 'o')) // ISMIP-HOM
        setExperName('1');
    else if ((temp == 'P') || (temp == 'p')) // ISMIP-POLICE
        setExperName('2');
  } else { // set a default: EISMINT II experiment A
    setExperName('A');
  }
  
  return 0;
}


PetscErrorCode IceEISModel::initFromOptions() {
  PetscErrorCode      ierr;
  char                inFile[PETSC_MAX_PATH_LEN];
  const PetscScalar   G_geothermal   = 0.042;             // J/m^2 s; geo. heat flux
  const PetscScalar   L              = 750e3;             // Horizontal extent of grid

  ierr = setExperNameFromOptions(); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  if (inFileSet == PETSC_TRUE) {
    ierr = initFromFile(inFile); CHKERRQ(ierr);
  } else { 
    ierr = initIceParam(grid.com, &grid.p, &grid.bag); CHKERRQ(ierr);
    grid.p->Mbz = 0;
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);
    // if no inFile then starts with zero ice
    ierr = VecSet(vh, 0);
    ierr = VecSet(vH, 0);
    ierr = VecSet(vbed, 0);
    ierr = VecSet(vGhf, G_geothermal);
    ierr = VecSet(vtau, DEFAULT_INITIAL_AGE_YEARS);
    // none use Goldsby-Kohlstedt
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    // all have no uplift at start
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);

    // next block checks if experiment is implemented
    // note height of grid must be great enough to handle max thickness
    switch (getExperName()) {
      case 'A':
        ierr = grid.rescale(L, L, 4500); CHKERRQ(ierr);
        break;
      case 'B':
      case 'C':
      case 'D':
        ierr = grid.rescale(L, L, 4000); CHKERRQ(ierr);
        break;
      case 'F':
        switch (getflowlawNumber()) {
          case 0:
          case 3:
            ierr = grid.rescale(L, L, 5000); CHKERRQ(ierr);
            break;
          case 1:
          case 4:
          case 5:
            ierr = grid.rescale(L, L, 6000); CHKERRQ(ierr);
            break;
          case 2:
            ierr = grid.rescale(L, L, 7000); CHKERRQ(ierr);
            break;
          default:  SETERRQ(1,"should not reach here (switch for rescale)\n");
        }
        break;
      case '0': // ISMIP-HEINO
        ierr = grid.rescale(2000e3, 2000e3, 5000); CHKERRQ(ierr);
        break;
      default:  SETERRQ(1,"ERROR: desired simplified geometry experiment NOT IMPLEMENTED\n");
    }
  }

  ierr = applyDefaultsForExperiment(); CHKERRQ(ierr);
  ierr = initAccumTs(); CHKERRQ(ierr);
  
  if (inFileSet == PETSC_FALSE) {
    ierr = fillintemps(); CHKERRQ(ierr);
  }
  
  ierr = afterInitHook(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::applyDefaultsForExperiment() {
  PetscErrorCode   ierr;
  PetscScalar      **mask;

  setThermalBedrock(PETSC_FALSE);
  setUseMacayealVelocity(PETSC_FALSE);
  setIsDrySimulation(PETSC_TRUE);
  setDoGrainSize(PETSC_FALSE);

  if (getExperName() == '0') { // ISMIP-HEINO
    setEnhancementFactor(3.0);
    setOceanKill(PETSC_TRUE);
    ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar myx = -grid.p->Lx + grid.p->dx * i, 
                          myy = -grid.p->Ly + grid.p->dy * j;
        const PetscScalar r = sqrt(PetscSqr(myx) + PetscSqr(myy));
        if (r > (2000e3 - 1))
          mask[i][j] = MASK_FLOATING_OCEAN0;
        else
          mask[i][j] = MASK_SHEET;
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  } else {
    setEnhancementFactor(1.0);
    ierr = VecSet(vMask, MASK_SHEET);
  }

  return 0;
}


PetscErrorCode IceEISModel::initAccumTs() {
  const PetscScalar   S_b     = 1e-2 * 1e-3 / secpera;    // Grad of accum rate change
  PetscScalar   S_T     = 1.67e-2 * 1e-3;           // K/m  Temp gradient
  PetscScalar         M_max,R_el,T_min;
  PetscErrorCode      ierr;
  PetscScalar         **accum, **Ts;

  switch (getExperName()) {
    case 'A':
      // start with zero ice and:
      M_max = 0.5 / secpera;  // Max accumulation
      R_el = 450e3;           // Distance to equil line (accum=0)
      T_min = 238.15;
      break;
    case 'B':
      // supposed to start from end of experiment A and:
      M_max = 0.5 / secpera;
      R_el = 450e3;
      T_min = 243.15;
      break;
    case 'C':
      // supposed to start from end of experiment A and:
      M_max = 0.25 / secpera;
      R_el = 425e3;
      T_min = 238.15;
      break;
    case 'D':
      // supposed to start from end of experiment A and:
      M_max = 0.5 / secpera;
      R_el = 425e3;
      T_min = 238.15;
      break;
    case 'F':
      // start with zero ice and:
      M_max = 0.5 / secpera;
      R_el = 450e3;
      T_min = 223.15;
      break;
    case '0':
      M_max = 0.3 / secpera;
      T_min = 233.15;
      S_T = 2.5e-9 / 1e9; // 2.5e-9 K km^-3
      break;
    default:
      SETERRQ(1,"\n experiment name unknown in IceEISModel::initAccumTs()\n");
  }
  globalMinTemp=T_min;

  // now fill in accum and surface temp
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar r = sqrt( PetscSqr(-grid.p->Lx + grid.p->dx*i)
                                  + PetscSqr(-grid.p->Ly + grid.p->dy*j) );
      // set accumulation
      if ((getExperName() >= 'A') && (getExperName() <= 'F'))
        accum[i][j] = PetscMin(M_max, S_b * (R_el-r));  // formula (7) in (Payne et al 2000)
      else if (getExperName() == '0') {
        const PetscScalar M_min = 0.15 / secpera;
        accum[i][j] = M_min + ((M_max - M_min) / 2000e3) * r;
      }
      // set surface temperature
      if ((getExperName() >= 'A') && (getExperName() <= 'F'))
        Ts[i][j] = T_min + S_T * r;                 // formula (8) in (Payne et al 2000)
      else if (getExperName() == '0')
        Ts[i][j] = T_min + S_T * r * r * r;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEISModel::fillintemps() {
  PetscErrorCode      ierr;
  PetscScalar         **Ts, ***T, ***Tb;

  // fill in all temps with Ts
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<grid.p->Mz; k++)
        T[i][j][k] = Ts[i][j];
      for (PetscInt k=0; k<grid.p->Mbz; k++)
        Tb[i][j][k] = Ts[i][j];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3b, vTb, INSERT_VALUES, vTb); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3b, vTb, INSERT_VALUES, vTb); CHKERRQ(ierr);
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid    g(com, rank, size);
    IceType*   ice;
    PetscInt   flowlawNumber = 0; // use Paterson-Budd by default
    
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);

    IceEISModel m(g, *ice);

    m.setflowlawNumber(flowlawNumber);
    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);
    ierr = m.setSoundingFromOptions(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "running simplified geometry experiment %c ...\n",m.getExperName()); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "done with run ... "); CHKERRQ(ierr);
    // see comments in run_ice.cc re default output naming convention
    ierr = m.writeFiles("simp_exper"); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
