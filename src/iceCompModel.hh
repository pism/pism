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

#include "materials.hh"
#include "iceModel.hh"


class IceCompModel : public IceModel {

public:
  IceCompModel(IceGrid &g, ThermoGlenArrIce &i);
  virtual ~IceCompModel();
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode initFromOptions();
  void setTest(char);
  PetscErrorCode setExactOnly(PetscTruth);
  PetscErrorCode summaryEismint_verify(bool);
  PetscErrorCode run();
  PetscErrorCode reportErrors();
  virtual PetscErrorCode dumpToFile_Matlab(const char *fname);

protected:
  virtual PetscErrorCode afterInitHook();
  PetscErrorCode createCompVecs();
  PetscErrorCode destroyCompVecs();
  PetscErrorCode createCompViewers();
  PetscErrorCode destroyCompViewers();
  PetscErrorCode updateCompViewers();

  // If I need features not present in IceType, I can get them here
  ThermoGlenArrIce &tgaIce;
  PetscTruth       testchosen, exactOnly, compVecsCreated, compViewersCreated;
  char             testname;  
  Vec              vSigmaComp;     // 3-D vector:   Mx x My x Mz
  PetscViewer      SigmaCompView, compSigmaMapView;

private:
  //general
  void mapcoords(const PetscInt i, const PetscInt j,
                 PetscScalar &x, PetscScalar &y, PetscScalar &r);
  virtual PetscScalar basalVelocity(const PetscScalar x, const PetscScalar y,
                                    const PetscScalar H, const PetscScalar T,
                                    const PetscScalar alpha, const PetscScalar mu);

  // tests A, B, C, D, E: separate public domain source "exactTestsABCDE.h/c"
  // test H: "exactTestH.h/c
  PetscErrorCode initTestISO();
  PetscErrorCode updateTestISO();

  // tests F & G: separate public domain source file "exactTestsFG.h/c" contains
  // "bothexact" which actually computes tests F & G 
  PetscErrorCode initTestFG();
  PetscErrorCode updateTestFG();

  PetscErrorCode computeGeometryErrors(         // alltests
        PetscScalar &gvolexact, PetscScalar &gareaexact, PetscScalar &gdomeHexact,
        PetscScalar &volerr, PetscScalar &areaerr,
        PetscScalar &gmaxHerr, PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
        PetscScalar &centerHerr);
  PetscErrorCode computeTemperatureErrors(      // tests F and G
        PetscScalar &gmaxTerr, PetscScalar &gavTerr);
  PetscErrorCode computeBasalTemperatureErrors( // tests F and G
        PetscScalar &gmaxTerr, PetscScalar &gavTerr, PetscScalar &centerTerr);
  PetscErrorCode computeSigmaErrors(            // tests F and G
        PetscScalar &gmaxSigmaerr, PetscScalar &gavSigmaerr);
  PetscErrorCode computeSurfaceVelocityErrors(  // tests F and G
        PetscScalar &gmaxUerr, PetscScalar &gavUerr,  // 2D vector errors
        PetscScalar &gmaxWerr, PetscScalar &gavWerr); // scalar errors
  PetscErrorCode computeBasalVelocityErrors(    // test E
        PetscScalar &exactmaxspeed,
        PetscScalar &gmaxvecerr, PetscScalar &gavvecerr,
        PetscScalar &gmaxuberr, PetscScalar &gmaxvberr);

private:
  // boundary conditions & parameters for all tests
  static PetscScalar Ggeo;    // J/m^2 s; geothermal heat flux, assumed constant
  static PetscScalar ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static PetscScalar Tmin;    // K;       minimum temperature (at center)
  static PetscScalar LforFG;  // m;  exact radius of tests F&G ice sheet
  static PetscScalar ApforG;  // m;  magnitude A_p of annular perturbation for test G;
                              //     note period t_p is set internally to 2000 years
  static PetscScalar ablationRateOutside;
  PetscScalar        f;       // ratio of ice density to bedrock density
};

#endif /* __iceCompModel_hh */
