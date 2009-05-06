// Copyright (C) 2004-2009 Nathan Shemonski and Ed Bueler
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

#ifndef __iceGRNModel_hh
#define __iceGRNModel_hh

#include <petscvec.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../coupler/pccoupler.hh"


//! Implements upper ice surface climate inputs, as parameterized in EISMINT-Greenland experiments.
/*!
Gives an elevation- and latitude-dependent mean annual surface temperature and 
amount of summer warming.

Also implements GWL3 climate warming scenario.
 */
class PISMEISGREENPDDCoupler : public PISMPDDCoupler {

public:
  PISMEISGREENPDDCoupler();
  PetscErrorCode startGWL3AtYear(PetscScalar year); // call with start year to do GML3
  virtual PetscErrorCode updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, IceModelVec2* &pvst);

protected:
  virtual PetscScalar getSummerWarming(
             const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma);
  virtual PetscScalar calculateMeanAnnual(PetscScalar h, PetscScalar lat);

private:
  PetscTruth  doGWL3;
  PetscScalar startGWL3Year;
};


//! Implements EISMINT-Greenland experiments.
/*!
This derived class adds, essentially, only the minimum functionality needed
to implement the choices stated in \ref RitzEISMINT , the EISMINT-Greenland 
specification:
- A PDD model is always used, through PISMEISGREENPDDCoupler.
- An enhancement factor of 3.0 is used.
- -ocean_kill is used by default.
- Constant geothermal flux.
- There is special code to ``clean out'' Ellsmere Island (and Iceland) so ice won't spread to
  edge of computational grid; this should probably be moved to the scripts which set up the
  bootstrap file.

A separate driver is used, namely src/pgrn.cc.
 */
class IceGRNModel : public IceModel {

public:
  IceGRNModel(IceGrid &g);
  PetscErrorCode attachEISGREENPDDPCC(PISMEISGREENPDDCoupler &p);
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode init_couplers();
  virtual PetscErrorCode set_vars_from_options();

protected:
  PISMEISGREENPDDCoupler *pddPCC; // points to same PISMAtmosCoupler as IceModel::atmosPCC,
                                  //   but we access PDD parameters through this pointer

  int expernum;  // SSL2 is 1, CCL3 is 3, GWL3 is 4

private:
  PetscTruth  noEllesmereIcelandDelete,
              haveSurfaceTemp,
              haveGeothermalFlux;
  PetscErrorCode ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat);
  PetscErrorCode cleanExtraLand();
};
#endif

