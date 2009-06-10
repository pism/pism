// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
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

#ifndef __iceEnthalpyModel_hh
#define __iceEnthalpyModel_hh

#include <petsc.h>
#include "iceModelVec.hh"
#include "materials.hh"
#include "columnSystem.hh"
#include "iceModel.hh"


//! Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
See \ref AschwandenBlatter2009.  The basic references are \ref PatersonBudd 
and \ref LliboutryDuval1985.
 */
class PolyThermalGlenPBLDIce : public ThermoGlenIce {
public:
  PolyThermalGlenPBLDIce(MPI_Comm c,const char pre[]);

  using ThermoGlenIce::flow;
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const;

  using ThermoGlenIce::effectiveViscosityColumn;
  virtual PetscScalar effectiveViscosityColumn(
                PetscScalar H,  PetscInt kbelowH, const PetscScalar *zlevels,
                PetscScalar u_x,  PetscScalar u_y, PetscScalar v_x,  PetscScalar v_y,
                const PetscScalar *T1, const PetscScalar *T2) const;

  /* these are used in src/base/ssaJed/ stuff only, so not re-implemented for now:
    integratedStoreSize(), integratedStore(), integratedViscosity()
  */

  NCConfigVariable *config;

protected:
  PetscReal water_frac_coeff;
};


// Tridiagonal linear system for vertical column of enthalpy-based conservation of energy.
/*
Call sequence like this:
\code
  enthSystemCtx foo;
  foo.dx = ...  // set public constants
  foo.u = ...   // set public pointers
  foo.initAllColumns();
  for (i in ownership) {
    for (j in ownership) {
      ks = ...
      foo.setIndicesThisColumn(i,j,ks);
      [COMPUTE OTHER PARAMS]
      foo.setSchemeParamsThisColumn(mask,isMarginal,lambda);  
      foo.setSurfaceBoundaryValuesThisColumn(Ts);
      foo.setBasalBoundaryValuesThisColumn(Ghf,Tshelfbase,Rb);
      foo.solveThisColumn(x);
    }  
  }
\endcode
 */
#if 0
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(int my_Mz, int my_Mbz);
  PetscErrorCode initAllColumns();
  PetscErrorCode setIndicesThisColumn(PetscInt i, PetscInt j, PetscInt ks);  
  PetscErrorCode setSchemeParamsThisColumn(
                     PetscScalar my_mask, bool my_isMarginal, PetscScalar my_lambda);  
  PetscErrorCode setSurfaceBoundaryValuesThisColumn(PetscScalar my_Ts);
  PetscErrorCode setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Tshelfbase, PetscScalar my_Rb);
  PetscErrorCode solveThisColumn(PetscScalar **x);  

public:
  // constants which should be set before calling initForAllColumns()
  PetscScalar  dx,
               dy,
               dtTemp,
               dzEQ,
               dzbEQ,
               ice_rho,
               ice_c_p,
               ice_k,
               bed_thermal_rho,
               bed_thermal_c_p,
               bed_thermal_k;
  // pointers which should be set before calling initForAllColumns()
  PetscScalar  *T,
               *Tb,
               *u,
               *v,
               *w,
               *Sigma;
  IceModelVec3 *T3;

protected: // used internally
  PetscInt    Mz, Mbz, k0;
  PetscInt    i, j, ks;
  PetscScalar mask, lambda, Ts, Ghf, Tshelfbase, Rb;
  bool        isMarginal;
  PetscScalar nuEQ,
              rho_c_I,
              rho_c_br,
              rho_c_av,
              iceK,
              iceR,
              brK,
              brR,
              rho_c_ratio,
              dzav,
              iceReff,
              brReff;
  bool        initAllDone,
              indicesValid,
              schemeParamsValid,
              surfBCsValid,
              basalBCsValid;
};
#endif


//! Temporary class for development of enthalpy-based polythermal PISM.
/*!
Based for now on Bueler's reading of \ref AschwandenBlatter2009.
 */
class IceEnthalpyModel : public IceModel {

public:
  IceEnthalpyModel(IceGrid &g);

  using IceModel::initFromFile;
  virtual PetscErrorCode initFromFile(const char *);

  using IceModel::write_extra_fields;
  virtual PetscErrorCode write_extra_fields(const char filename[]);

  bool doColdIceTemperatureStep;

protected:
  using IceModel::createVecs;
  virtual PetscErrorCode createVecs();
  
  using IceModel::init_physics;
  virtual PetscErrorCode init_physics();

  // PetscErrorCode setEnth3toCTSValue();
  virtual PetscErrorCode setEnth3FromT3_ColdIce();

  using IceModel::temperatureStep;
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);

  virtual PetscErrorCode enthalpyStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  
  using IceModel::temperatureAgeStep;
  virtual PetscErrorCode temperatureAgeStep();

  IceModelVec3  Enth3, EnthNew3;
};

#endif

