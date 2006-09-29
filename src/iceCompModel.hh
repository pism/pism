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
  PetscErrorCode reporterror();
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
  PetscTruth       testchosen;
  char             testname;  
  PetscTruth       exactOnly;
  Vec              vSigmaComp;     // 3-D vector:   Mx x My x Mz
  PetscViewer      SigmaCompView, compSigmaMapView;

private:
  //general
  void mapcoords(const PetscInt i, const PetscInt j,
                 PetscScalar &x, PetscScalar &y, PetscScalar &r);

  // tests B, C, D & H: separate public domain source "exactTestsBCD.cc" contains "exactB", 
  // "exactC", "exactD" to compute these tests, along with "exactH" here.
  PetscErrorCode initTestBCDH();
  PetscErrorCode updateTestBCDH();

  // Test H is the simple isostasy version of test C followed by a version of test B at
  // t0=40034 years.  It is used in
  //    Bueler and others (2006) "Fast computation of a viscoelastic deformable Earth model
  //    for ice sheet simulations", submitted to Ann. Glaciol.
  int exactH(double t, double r, double &H, double &M);
  
  // tests F & G: separate public domain source file "exactTestsFG.cc" contains "bothexact"
  // which actually computes tests F & G 
  PetscErrorCode initTestFG();
  PetscErrorCode updateTestFG();

private:
  // boundary condition parameters for all tests
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
