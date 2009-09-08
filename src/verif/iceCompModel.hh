// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceCompModel_hh
#define __iceCompModel_hh

#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../base/iceModelVec.hh"

class IceCompModel : public IceModel {

public:
  IceCompModel(IceGrid &g, int mytest);
  virtual ~IceCompModel();
  
  // re-defined steps of init() sequence:
  virtual PetscErrorCode set_grid_defaults();     // called by IceModel::grid_setup()
  virtual PetscErrorCode set_grid_from_options(); // called by IceModel::grid_setup()
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode init_physics();
  virtual PetscErrorCode init_couplers();
  virtual PetscErrorCode set_vars_from_options(); // called by IceModel::model_state_setup()

  virtual PetscErrorCode createViewers();
  virtual PetscErrorCode destroyViewers();

  PetscErrorCode reportErrors();

protected:
  // related to all (or most) tests
  ThermoGlenArrIce *tgaIce;
  PetscTruth   exactOnly;
  int          testname;
  void         mapcoords(const PetscInt i, const PetscInt j,
                  PetscScalar &x, PetscScalar &y, PetscScalar &r);
  virtual PetscErrorCode 
               summaryPrintLine(PetscTruth printPrototype, bool tempAndAge,
                   PetscScalar year,  PetscScalar dt, 
                   PetscScalar volume_kmcube,  PetscScalar area_kmsquare,
                   PetscScalar meltfrac,  PetscScalar H0,  PetscScalar T0);
  virtual PetscErrorCode 
               additionalAtStartTimestep();
  virtual PetscErrorCode 
               additionalAtEndTimestep();
  PetscErrorCode computeGeometryErrors(    // all tests except K
        PetscScalar &gvolexact, PetscScalar &gareaexact, PetscScalar &gdomeHexact,
        PetscScalar &volerr, PetscScalar &areaerr,
        PetscScalar &gmaxHerr, PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
        PetscScalar &centerHerr);

  // related to tests A B C D E H
  PetscErrorCode initTestABCDEH();
  PetscErrorCode getCompSourcesTestCDH();  // only for time-dependent compensatory sources
  PetscErrorCode fillSolnTestABCDH();  // only used with exactOnly == PETSC_TRUE
  
  // related to test E
  PetscErrorCode fillSolnTestE();  // only used with exactOnly == PETSC_TRUE
  PetscErrorCode computeBasalVelocityErrors(    // test E only
        PetscScalar &exactmaxspeed,
        PetscScalar &gmaxvecerr, PetscScalar &gavvecerr,
        PetscScalar &gmaxuberr, PetscScalar &gmaxvberr);
  virtual PetscScalar basalVelocitySIA( // not recommended, generally
                             PetscScalar x, PetscScalar y, PetscScalar H, PetscScalar T,
                             PetscScalar alpha, PetscScalar mu, PetscScalar min_T) const;

  // related to test L
  IceModelVec2   vHexactL;
  PetscErrorCode initTestL();
  PetscErrorCode fillSolnTestL();  // only used with exactOnly == PETSC_TRUE

  // related to tests F G; see iCMthermo.cc
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  PetscErrorCode initTestFG();
  PetscErrorCode getCompSourcesTestFG();
  PetscErrorCode fillSolnTestFG();  // only used with exactOnly == PETSC_TRUE
  PetscErrorCode computeTemperatureErrors(      // tests F and G
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr);
  PetscErrorCode computeBasalTemperatureErrors( // tests F and G
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr, PetscScalar &centerTerr);
  PetscErrorCode computeSigmaErrors(            // tests F and G
                   PetscScalar &gmaxSigmaerr, PetscScalar &gavSigmaerr);
  PetscErrorCode computeSurfaceVelocityErrors(  // tests F and G
                   PetscScalar &gmaxUerr, PetscScalar &gavUerr,  // 2D vector errors
                   PetscScalar &gmaxWerr, PetscScalar &gavWerr); // scalar errors
  
  IceModelVec3   SigmaComp3;

  PetscViewer    SigmaCompView, compSigmaMapView;
  PetscErrorCode updateCompViewers();

  // related to test K; see iCMthermo.cc
  PetscErrorCode initTestK();
  PetscErrorCode fillSolnTestK();  // only used with exactOnly == PETSC_TRUE
  PetscErrorCode computeIceBedrockTemperatureErrors( // test K only
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr,
                   PetscScalar &gmaxTberr, PetscScalar &gavTberr);


private:
  PetscScalar        f;       // ratio of ice density to bedrock density
  PetscTruth         bedrock_is_ice_forK;

  static const PetscScalar ablationRateOutside;

  // see iCMthermo.cc
  static const PetscScalar Ggeo;    // J/m^2 s; geothermal heat flux, assumed constant
  static const PetscScalar ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static const PetscScalar Tmin;    // K;       minimum temperature (at center)
  static const PetscScalar LforFG;  // m;  exact radius of tests F&G ice sheet
  static const PetscScalar ApforG;  // m;  magnitude A_p of annular perturbation for test G;
  // period t_p is set internally to 2000 years
};

#endif /* __iceCompModel_hh */
