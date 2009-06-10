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
#include "iceModel.hh"
#include "iceModelVec.hh"


//! Temporary class for development of enthalpy-based polythermal PISM.
/*!
Based for now on Bueler's reading of A. Aschwandedn and H. Blatter,
"An enthalpy formulation for polythermal glaciers and ice sheets", in preparation.
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
  
  // PetscErrorCode setEnth3toCTSValue();
  virtual PetscErrorCode setEnth3FromTemp_ColdIce();

  using IceModel::temperatureStep;
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);

  virtual PetscErrorCode enthalpyStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  
  using IceModel::temperatureAgeStep;
  virtual PetscErrorCode temperatureAgeStep();

  IceModelVec3  Enth3, EnthNew3;

private:

  PetscScalar getPressureFromDepth(PetscScalar depth);

  //! Get pressure melting temperature (K) and enthalpy at phase transitions from pressure p.
  /*! From \ref AschwandenBlatter2009,
     \f[ T_m(p) = T_0 - \beta p, \f]
     \f[ H_l(p) = c_w T_m(p), \f]
     \f[ H_s(p) = -L + H_l(p). \f]
  Also returns H_s.
   */ 
  PetscScalar get_H_s(PetscScalar p, PetscScalar &T_m, PetscScalar &H_l, PetscScalar &H_s);

  //! Get absolute ice temperature (K) from enthalpy H and pressure p.
  /*! From \ref AschwandenBlatter2009,
     \f[ T=T(H,p) = \begin{cases} c_i^{-1} (H-H_s(p)) + T_m(p), & H < H_s(p), \\
                                  T_m(p), &                       H_s(p) \le H < H_l(p). \end{cases} \f]

  We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
  \f$H \ge H_l(p)\f$.
  
  See get_H_s() for calculation of \f$T_m(p)\f$, \f$H_l(p)\f$, and \f$H_s(p)\f$. 
   */ 
  PetscScalar getAbsTemp(PetscScalar H, PetscScalar p);

  //! Get liquid water fraction from enthalpy H and pressure p.
  /*! From \ref AschwandenBlatter2009,
     \f[ \omega=\omega(H,p) = \begin{cases} 0.0,            & H \le H_s(p), \\
                                            (H-H_s(p)) / L, & H_s(p) < H < H_l(p). \end{cases} \f]

  We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
  \f$H \ge H_l(p)\f$.
  
  See get_H_s() for calculation of \f$T_m(p)\f$, \f$H_l(p)\f$, and \f$H_s(p)\f$. 
   */ 
  PetscScalar getWaterFraction(PetscScalar H, PetscScalar p);

  PetscScalar getEnth(PetscScalar T, PetscScalar omega, PetscScalar p);

  PetscScalar getEnth_pa(PetscScalar T_pa, PetscScalar omega, PetscScalar p);

};

#endif

