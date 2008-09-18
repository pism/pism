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

#ifndef __iceCompModel_hh
#define __iceCompModel_hh

#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../base/iceModelVec.hh"

class IceCompModel : public IceModel {

public:
  IceCompModel(IceGrid &g, IceType *i, const char mytest);
  virtual ~IceCompModel();
  using IceModel::initFromOptions;
  virtual PetscErrorCode initFromOptions();
  PetscErrorCode reportErrors();

protected:
  // related to all (or most) tests
  ThermoGlenArrIce *tgaIce;
  PetscTruth   exactOnly;
  char         testname;
  void         mapcoords(const PetscInt i, const PetscInt j,
                  PetscScalar &x, PetscScalar &y, PetscScalar &r);
  virtual PetscErrorCode 
               summaryPrintLine(const PetscTruth printPrototype, const PetscTruth tempAndAge,
                  const PetscScalar year, const PetscScalar dt, 
                  const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
                  const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);
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
  virtual PetscScalar basalVelocity(const PetscScalar x, const PetscScalar y,
                                    const PetscScalar H, const PetscScalar T,
                                    const PetscScalar alpha, const PetscScalar mu);

  // related to test L
  Vec            vHexactL;
  PetscTruth     vHexactLCreated;
  PetscErrorCode initTestL();
  PetscErrorCode fillSolnTestL();  // only used with exactOnly == PETSC_TRUE

  // related to tests F G; see iCMthermo.cc
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount);
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
  PetscTruth     compVecsCreated, compViewersCreated;
  PetscErrorCode createCompVecs();
  PetscErrorCode destroyCompVecs();
  PetscErrorCode createCompViewers();
  PetscErrorCode destroyCompViewers();
  PetscErrorCode updateCompViewers();

  // related to test K; see iCMthermo.cc
  PetscErrorCode initTestK();
  PetscErrorCode fillSolnTestK();  // only used with exactOnly == PETSC_TRUE
  PetscErrorCode computeIceBedrockTemperatureErrors( // test K only
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr,
                   PetscScalar &gmaxTberr, PetscScalar &gavTberr);


private:
  static PetscScalar ablationRateOutside;
  PetscScalar        f;       // ratio of ice density to bedrock density
  PetscTruth         bedrock_is_ice_forK;

  // see iCMthermo.cc
  static PetscScalar Ggeo;    // J/m^2 s; geothermal heat flux, assumed constant
  static PetscScalar ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static PetscScalar Tmin;    // K;       minimum temperature (at center)
  static PetscScalar LforFG;  // m;  exact radius of tests F&G ice sheet
  static PetscScalar ApforG;  // m;  magnitude A_p of annular perturbation for test G;
                              //     note period t_p is set internally to 2000 years
};

#endif /* __iceCompModel_hh */
