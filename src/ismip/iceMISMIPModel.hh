// Copyright (C) 2008 Ed Bueler
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

#ifndef __iceMISMIPModel_hh
#define __iceMISMIPModel_hh

#include <petscvec.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"

class MISMIPIce : public ThermoGlenIce {
public:
  MISMIPIce();
  using ThermoGlenIce::flow;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  // this one returns nu; ignors temp & pressure
  virtual PetscScalar effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const;
  // this one returns nu * H; calls effectiveViscosity() for nu
  virtual PetscScalar effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscInt kbelowH,
                           const PetscInt nlevels, PetscScalar *zlevels,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const;
  PetscErrorCode setA(const PetscScalar myA);
  PetscScalar    softnessParameter(const PetscScalar T) const;
  PetscScalar    hardnessParameter(const PetscScalar T) const;
  PetscErrorCode printInfo(const int thresh, MPI_Comm com);
  
protected:
  PetscScalar  A_MISMIP, // softness
               B_MISMIP; // hardness; B = A^{-1/n}
};


struct routineStatsType {
  PetscScalar jg, xg, hxg, maxubar, avubarG, avubarF, dHdtnorm;
};


struct mismipStatsType {
  PetscScalar dxgdt, x1, x2, x3, h1, h2, h3, b1, b2, b3, q1, q2, q3;
};


//! Derived class of IceModel with performs MISMIP experiments.
/*!
See \e User's \e Manual and run script   examples/mismip/mismip.sh.

WE'VE GOT A PROBLEM WITH FLUX ACROSS GROUNDING LINE; not surprising ...
See "cflx" in output file.

I think the underlying issue *may* still be *regularity of the basal shear stress*, 
that is, of the shear stress coefficient, which jumps from C_MISMIP to zero across
the grounding line.

I did an experiment that made the ice grounded
all the way out to the calving front, but has zero basal resistance; the results are
similar the intended MISMIP case, even though the surface is not at all what is given
by the floatation criterion.

Damping out C_MISMIP in the 50km on the grounded side of the grounding line (a linear
decrease linearly from C_MISMIP at xg - 50km to zero at xg) does not make a difference,
really.  In other words, I thought we needed
   \tau_b \in W^{1,\infty} 
or something; without the damping all we have is
   \tau_b \in L^\infty.
The results in the former case, of smoother C_MISMIP, are no better, I think.
 */
class IceMISMIPModel : public IceModel {

public:
  IceMISMIPModel(IceGrid &g, MISMIPIce *mismip_i);
  virtual ~IceMISMIPModel(); // must be virtual merely because some members are virtual

  virtual PetscErrorCode setFromOptions();
  using IceModel::initFromOptions;
  virtual PetscErrorCode initFromOptions();
  PetscErrorCode         additionalAtStartTimestep();
  PetscErrorCode         additionalAtEndTimestep();
  virtual PetscErrorCode summaryPrintLine(
                const PetscTruth printPrototype, const PetscTruth tempAndAge,
                const PetscScalar year, const PetscScalar dt, 
                const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
                const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);
  virtual PetscScalar    basalDragx(PetscScalar **tauc,
                                    PetscScalar **u, PetscScalar **v,
                                    PetscInt i, PetscInt j) const;
  virtual PetscScalar    basalDragy(PetscScalar **tauc,
                                    PetscScalar **u, PetscScalar **v,
                                    PetscInt i, PetscInt j) const;
  PetscErrorCode         printBasalAndIceInfo();
  PetscErrorCode         writeMISMIPFinalFiles();

private:
  MISMIPIce   *mismip_ice;
  PetscInt    exper, gridmode, stepindex, modelnum;
  char        sliding;
  PetscScalar runtimeyears, dHdtnorm_atol;
  PetscTruth  writeExtras, steadyOrGoalAchieved;
  char        initials[PETSC_MAX_PATH_LEN],  // initials of user, for MISMIP reporting
              mprefix[PETSC_MAX_PATH_LEN];
  PetscViewer tviewfile;
  char        tfilename[PETSC_MAX_PATH_LEN];

  routineStatsType  rstats;
  mismipStatsType   mstats;

  PetscErrorCode  setMISMIPBed();
  PetscErrorCode  setMISMIPMask();
  PetscErrorCode  getMISMIPStats();
  PetscErrorCode  getRoutineStats();

  PetscScalar m_MISMIP, C_MISMIP;
  PetscScalar regularize_MISMIP;

  PetscScalar basalIsotropicDrag(PetscScalar **u, PetscScalar **v, 
                                 PetscInt i, PetscInt j) const;
  PetscErrorCode writeMISMIPasciiFile(const char mismiptype, char* filename);
};

#endif  // __iceMISMIPModel_hh

