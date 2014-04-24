// Copyright (C) 2009-2011, 2013, 2014 Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef __tempSystem_hh
#define __tempSystem_hh

#include <petscsys.h>

#include "columnSystem.hh"
#include "pism_const.hh"
#include "Mask.hh"

namespace pism {

class IceModelVec3;

//! Tridiagonal linear system for vertical column of temperature-based conservation of energy problem.
/*!
  Call sequence like this:
  \code
  tempSystemCtx foo;
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
class tempSystemCtx : public columnSystemCtx {

public:
  tempSystemCtx(int my_Mz, std::string my_prefix);
  PetscErrorCode initAllColumns();
  PetscErrorCode setSchemeParamsThisColumn(MaskValue my_mask, bool my_isMarginal, double my_lambda);  
  PetscErrorCode setSurfaceBoundaryValuesThisColumn(double my_Ts);
  PetscErrorCode setBasalBoundaryValuesThisColumn(
                                                  double my_G0, double my_Tshelfbase, double my_Rb);

  PetscErrorCode solveThisColumn(double *x);  

public:
  // constants which should be set before calling initForAllColumns()
  double  dx,
    dy,
    dtTemp,
    dzEQ,
    ice_rho,
    ice_c_p,
    ice_k;
  // pointers which should be set before calling initForAllColumns()
  double  *T,
    *u,
    *v,
    *w,
    *strain_heating;
  IceModelVec3 *T3;

protected: // used internally
  int    Mz;
  double lambda, Ts, G0, Tshelfbase, Rb;
  MaskValue    mask;
  bool        isMarginal;
  double nuEQ,
    rho_c_I,
    iceK,
    iceR;
  bool        initAllDone,
    indicesValid,
    schemeParamsValid,
    surfBCsValid,
    basalBCsValid;
};

} // end of namespace pism

#endif  /* __tempSystem_hh */

