// Copyright (C) 2004-2009 Nathan Shemonski, Ed Bueler and Constantine Khroulev
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
#include "../base/iceModel.hh"

//! Very slightly modified mass balance formulas for EISMINT-Greenland.
/*!
EISMINT-Greenland [\ref RitzEISMINT] has a slightly different interpretation
of positive degree day factors, compared to [\ref Faustoetal2009].  In particular,
we stop using formula (6) in [\ref Faustoetal2009].
 */
class EISGREENMassBalance : public PDDMassBalance {

public:
  EISGREENMassBalance(const NCConfigVariable& myconfig);

  //! Leaves pddFactorIce and pddFactorSnow alone; different from base class action.
  virtual PetscErrorCode setDegreeDayFactorsFromSpecialInfo(
             PetscScalar latitude, PetscScalar T_mj);
};


//! Upper ice surface climate inputs, as parameterized in EISMINT-Greenland experiments.
/*!
Gives an elevation- and latitude-dependent mean annual ice surface temperature and 
snow temperature, for PDD.  Also implements a greenhouse climate warming scenario,
for experiment GWL3 [\ref RitzEISMINT].
 */
class EISGREENAtmosCoupler : public PISMSnowModelAtmosCoupler {

public:
  EISGREENAtmosCoupler();

  using PISMSnowModelAtmosCoupler::initFromOptions; // overrides but calls this
  //! Just checks correct attachment of new mass balance scheme.
  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);

  using PISMSnowModelAtmosCoupler::updateSurfTempAndProvide; // overrides this
  //! Implements ice surface temperature parameterization.
  virtual PetscErrorCode updateSurfTempAndProvide(
               PetscScalar t_years, PetscScalar dt_years,
               IceModelVec2* &pvst);

  //! Turns on and gives starting model year for greenhouse warming scenario in experiment GWL3.
  PetscErrorCode startGreenhouseAtYear(PetscScalar year);

protected:

  using PISMSnowModelAtmosCoupler::parameterizedUpdateSnowSurfaceTemp; // overrides this
  //! Implements snow surface temperature parameterization, for use in mass balance model.
  virtual PetscErrorCode parameterizedUpdateSnowSurfaceTemp(
            PetscScalar t_years, PetscScalar dt_years);

  //! Implements mean annual temperature parameterization, used by parameterizedUpdateSnowSurfaceTemp().
  PetscScalar meanAnnualTemp(PetscScalar h, PetscScalar lat);

  //! When the greenhouse climate scenario is on this computes the scalar warming amount.
  PetscScalar shiftForGreenhouse(PetscScalar t_years, PetscScalar dt_years);

private:
  PetscTruth  doGreenhouse;
  PetscScalar startYearGreenhouse;
};


typedef enum {SSL2, SSL3, CCL3, GWL3} EISGREENrun;


//! Implements EISMINT-Greenland experiments.
/*!
This derived class adds only the minimum functionality needed
to implement the choices stated in \ref RitzEISMINT , the EISMINT-Greenland 
specification:
- A PDD model is always used, through a new coupler (EISGREENAtmosCoupler)
  calling a new mass balance scheme (EISGREENMassBalance).  These are small modifications
  of the default scheme from \ref Faustoetal2009, implemented in PISMSnowModelAtmosCoupler.
- An enhancement factor of 3.0 is used.
- -ocean_kill is used by default.

Which experiment to do is chosen by one of options -ssl2,-ccl3,-gwl3.  (Experiment
SSL3 is not implemented; see User's Manual.)

A separate driver calls this derived class and the new coupler, namely src/pgrn.cc.
 */
class IceGRNModel : public IceModel {

public:
  IceGRNModel(IceGrid &g) : IceModel(g) {}
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode init_couplers();

private:
  EISGREENrun exper;
};
#endif

