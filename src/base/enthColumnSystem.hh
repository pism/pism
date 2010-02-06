// Copyright (C) 2009-2010 Andreas Aschwanden and Ed Bueler
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

#ifndef __enthColumnSystem_hh
#define __enthColumnSystem_hh

#include "columnSystem.hh"

//! Tridiagonal linear system for vertical column of enthalpy-based conservation of energy.
/*!
See the page documenting \ref bombproofenth.  We implement equations
 FIXME: 
(\ref bombtwo), (\ref bedrockeqn),
(\ref geothermalbedeqn), (\ref icebedfinalcold), (\ref icebedfinaltemperate),
(\ref neartopofbedrock), and (\ref icebasenobedrock).

There are Mbz bedrock
temperature variables, a heat flux variable X at ice/bed interface,
and Mz ice enthalpy variables for a total of Mbz+Mz+1 variables.  These are
the variables at the corresponding levels:
  - x[0] .. x[Mbz-2]         temperature in bedrock
  - x[Mbz-1]                 temp at top of the bedrock layer; location z=0;
                             Mbz==1 \e is allowed, so may equal x[0]
  - x[Mbz]                   = X, a heat flux at the ice bed interface, with
                             different interpretations in different cases; see
                             tables in \ref bombproofenth
  - x[Mbz+1]                 enthalpy at base of ice; z=0 so same \e location
                             as [Mbz-1]
  - x[Mbz+1+1]..x[Mbz+1+ks]  enthalpy in the ice
  - x[Mbz+1+ks+1]..x[Mbz+Mz] enthalpy in the air

The equations have different meanings and different units at the different
levels.  Let "CE" stand for conservation of energy
  - eqn 0                    heat flux BC; units of Kelvin if at bottom of
                             but units of heat flux if Mbz==1
  - eqns 1 .. Mbz-2          CE in bedrock; units of Kelvin
  - eqn Mbz-1                heat flux BC; units of heat flux; could be eqn 0
  - eqn Mbz                  continuity of temperature at ice-bed interface;
                             units of Kelvin
  - eqn Mbz+1                heat flux BC; units of heat flux
  - eqns Mbz+1+1..Mbz+1+ks   CE in ice; units of enthalpy (J kg-1)
  - eqns Mbz+1+ks+1..Mbz+Mz  sets air enthalpy; units of enthalpy (J kg-1)

\e All levels are solved.  That is, the tridiagonal linear system \e always
has Mbz+Mz+1 equations, and the enthalpy of the air above the ice is always
solved-for.

Coefficients in the system (D,L,U expressions) are unitless.
*/
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(int my_Mz, int my_Mbz);

  PetscErrorCode initAllColumns();
  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);

  PetscErrorCode setSchemeParamsThisColumn(
            const bool my_isfloating, const bool my_ismarginal,
            const PetscScalar my_lambda);  
  PetscErrorCode setBoundaryValuesThisColumn(
            const PetscScalar my_Enth_surface,
            const PetscScalar my_Ghf, const PetscScalar my_Rb);

  PetscErrorCode solveThisColumn(PetscScalar **x);

public:
  // constants which should be set before calling initAllColumns()
  PetscScalar  dx,
               dy,
               dtTemp,
               dzEQ,
               dzbEQ,
               ice_rho,
               ice_c,
               ice_k,
               ice_nu,
               bed_rho,
               bed_c,
               bed_k;
  // pointers which should be set before calling initAllColumns()
  PetscScalar  *Enth,   // enthalpy in ice
               *Enth_s, // enthalpy level for CTS; function only of pressure
               *Tb,     // temperature in bedrock
               *u,
               *v,
               *w,
               *Sigma;
  IceModelVec3 *Enth3;

private: // used internally
  PetscInt    Mz, Mbz;
  PetscScalar lambda, nuEQ,
              iceK, iceRcold, iceRtemp,
              bedK, bedR,
              Enth_ks, Ghf, Rb;
  bool        initAllDone, schemeParamsValid, BCsValid,
              isfloating, ismarginal;
};

#endif

