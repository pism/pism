// Copyright (C) 2009-2011, 2013 Andreas Aschwanden and Ed Bueler
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

#ifndef __enthSystem_hh
#define __enthSystem_hh

#include <vector>

#include "columnSystem.hh"

class NCConfigVariable;
class IceModelVec3;
class EnthalpyConverter;

//! Tridiagonal linear system for conservation of energy in vertical column of ice enthalpy.
/*!
  See the page documenting \ref bombproofenth.  The top of
  the ice has a Dirichlet condition.
*/
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(const NCConfigVariable &config,
                IceModelVec3 &my_Enth3,
                PetscScalar my_dx,  PetscScalar my_dy,
                PetscScalar my_dt,  PetscScalar my_dz,
                int my_Mz, std::string my_prefix,
                EnthalpyConverter *my_EC);
  virtual ~enthSystemCtx();

  PetscErrorCode initThisColumn(int i, int j, bool my_ismarginal,
                                PetscReal ice_thickness,
                                PetscReal till_water_thickness,
                                IceModelVec3 *u3,
                                IceModelVec3 *v3,
                                IceModelVec3 *w3,
                                IceModelVec3 *strain_heating3);

  PetscScalar k_from_T(PetscScalar T);

  PetscErrorCode setDirichletSurface(PetscScalar my_Enth_surface);
  PetscErrorCode setDirichletBasal(PetscScalar Y);
  PetscErrorCode setBasalHeatFlux(PetscScalar hf);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);
  PetscErrorCode viewSystem(PetscViewer viewer) const;

  PetscErrorCode solveThisColumn(PetscScalar **x);

  int ks()
  { return m_ks; }

  PetscReal lambda()
  { return m_lambda; }
public:
  // arrays must be filled before calling solveThisColumn():
  PetscScalar  *Enth,   // enthalpy in ice at prev time step
    *Enth_s; // enthalpy level for CTS; function only of pressure
protected:
  PetscScalar *u, *v, *w, *strain_heating;

  PetscInt Mz, m_ks;

  PetscReal ice_rho, ice_c, ice_k, ice_K, ice_K0, p_air,
    dx, dy, dz, dt, nu, R_cold, R_temp, R_factor;

  PetscReal ice_thickness,
    m_lambda,              //!< implicit FD method parameter
    Enth_ks;             //!< top surface Dirichlet B.C.
  PetscReal a0, a1, b;   // coefficients of the first (basal) equation
  bool ismarginal, c_depends_on_T, k_depends_on_T;
  std::vector<PetscScalar> R; // values of k \Delta t / (\rho c \Delta x^2)

  IceModelVec3 *Enth3;
  EnthalpyConverter *EC;  // conductivity has known dependence on T, not enthalpy

  PetscErrorCode compute_enthalpy_CTS();
  PetscReal compute_lambda();

  virtual PetscErrorCode assemble_R();
  PetscErrorCode checkReadyToSolve();
};

#endif   //  ifndef __enthSystem_hh
