// Copyright (C) 2008-2010 Ed Bueler and Constantine Khroulev
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


struct routineStatsType {
  PetscScalar ig,		// index corresponding to the grounding line
    xg,				// grounding line position
    hxg,			// ice thickness at the grounding line
    maxubar,			// maximum vertically-averaged horizontal velocity
    avubarG,			// average of ubar over grounded points
    avubarF, 			// ditto, floating points
    dHdtnorm,			// maximum of the thickness rate of change 
    Ngrounded,			// number of grounded points
    Nfloating;			// number of floating points
};


struct mismipStatsType {
  PetscScalar dxgdt,		// rate of change of the grounding line position
    x1, x2, x3,			// grounding line position and its neighboring points
    h1, h2, h3,			// ice thicknesses corresponding to x1, x2, x3
    b1, b2, b3,			// bed depth below sea level at x1, x2, x3
    q1, q2, q3;			// signed ice flux at x1, x2, x3
};


//! Derived class of IceModel which performs MISMIP experiments.
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
\f[   \tau_b \in W^{1,\infty} \f]
or something; without the damping all we have is
\f[   \tau_b \in L^\infty. \f]
The results in the former case, of smoother C_MISMIP, are no better, I think.
 */
class IceMISMIPModel : public IceModel {

public:
  IceMISMIPModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides);
  virtual ~IceMISMIPModel(); // must be virtual merely because some members are virtual

  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode set_grid_defaults();
  virtual PetscErrorCode set_grid_from_options();
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode init_physics();
  virtual PetscErrorCode init_couplers();
  virtual PetscErrorCode set_time_from_options();
  virtual PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode misc_setup();
  

  PetscErrorCode         additionalAtStartTimestep();
  PetscErrorCode         additionalAtEndTimestep();
  virtual PetscErrorCode summaryPrintLine(
                 PetscTruth printPrototype,  bool tempAndAge,
                 PetscScalar year,  PetscScalar dt, 
                 PetscScalar volume,  PetscScalar area,
                 PetscScalar meltfrac,  PetscScalar H0,  PetscScalar T0);
  using IceModel::basalDragx;
  virtual PetscScalar    basalDragx(PetscScalar **tauc,
                                    PetscScalar **u, PetscScalar **v,
                                    PetscInt i, PetscInt j) const;
  using IceModel::basalDragy;
  virtual PetscScalar    basalDragy(PetscScalar **tauc,
                                    PetscScalar **u, PetscScalar **v,
                                    PetscInt i, PetscInt j) const;
  PetscErrorCode         printBasalAndIceInfo();
  PetscErrorCode         writeMISMIPFinalFiles();

private:
  PetscInt    exper, gridmode, stepindex, modelnum;
  char        sliding;
  PetscScalar initialthickness, runtimeyears, dHdtnorm_atol;
  PetscTruth  steadyOrGoalAchieved;
  bool writeExtras, tryCalving;
  char        initials[PETSC_MAX_PATH_LEN],  // initials of user, for MISMIP reporting
              mprefix[PETSC_MAX_PATH_LEN];
  PetscViewer tviewfile;
  char        tfilename[PETSC_MAX_PATH_LEN];

  routineStatsType  rstats;
  mismipStatsType   mstats;

  PetscErrorCode  setBed();
  PetscErrorCode  setMask();
  PetscErrorCode  getMISMIPStats();
  PetscErrorCode  getRoutineStats();

  PetscScalar m_MISMIP, C_MISMIP;
  PetscScalar regularize_MISMIP;

  PetscErrorCode calving();
  PetscScalar    basalIsotropicDrag(PetscScalar **u, PetscScalar **v, 
                                    PetscInt i, PetscInt j) const;
  PetscErrorCode writeMISMIPasciiFile(const char mismiptype, char* filename);
};

#endif  // __iceMISMIPModel_hh

