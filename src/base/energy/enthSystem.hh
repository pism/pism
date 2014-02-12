// Copyright (C) 2009-2011, 2013, 2014 Andreas Aschwanden and Ed Bueler
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

class PISMConfig;
class IceModelVec3;
class EnthalpyConverter;

//! Tridiagonal linear system for conservation of energy in vertical column of ice enthalpy.
/*!
  See the page documenting \ref bombproofenth.  The top of
  the ice has a Dirichlet condition.
*/
class enthSystemCtx : public columnSystemCtx {

public:
  enthSystemCtx(const PISMConfig &config,
                IceModelVec3 &my_Enth3,
                double my_dx,  double my_dy,
                double my_dt,  double my_dz,
                int my_Mz, std::string my_prefix,
                EnthalpyConverter *my_EC);
  virtual ~enthSystemCtx();

  PetscErrorCode initThisColumn(int i, int j, bool my_ismarginal,
                                double ice_thickness,
                                double till_water_thickness,
                                IceModelVec3 *u3,
                                IceModelVec3 *v3,
                                IceModelVec3 *w3,
                                IceModelVec3 *strain_heating3);

  double k_from_T(double T);

  PetscErrorCode setDirichletSurface(double my_Enth_surface);
  PetscErrorCode setDirichletBasal(double Y);
  PetscErrorCode setBasalHeatFlux(double hf);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);
  PetscErrorCode viewSystem(PetscViewer viewer) const;

  PetscErrorCode solveThisColumn(double *x);

  int ks()
  { return m_ks; }

  double lambda()
  { return m_lambda; }
public:
  // arrays must be filled before calling solveThisColumn():
  double  *Enth,   // enthalpy in ice at prev time step
    *Enth_s; // enthalpy level for CTS; function only of pressure
protected:
  double *u, *v, *w, *strain_heating;

  int Mz, m_ks;

  double ice_rho, ice_c, ice_k, ice_K, ice_K0, p_air,
    dx, dy, dz, dt, nu, R_cold, R_temp, R_factor;

  double ice_thickness,
    m_lambda,              //!< implicit FD method parameter
    Enth_ks;             //!< top surface Dirichlet B.C.
  double a0, a1, b;   // coefficients of the first (basal) equation
  bool ismarginal, c_depends_on_T, k_depends_on_T;
  std::vector<double> R; // values of k \Delta t / (\rho c \Delta x^2)

  IceModelVec3 *Enth3;
  EnthalpyConverter *EC;  // conductivity has known dependence on T, not enthalpy

  PetscErrorCode compute_enthalpy_CTS();
  double compute_lambda();

  virtual PetscErrorCode assemble_R();
  PetscErrorCode checkReadyToSolve();
};

#endif   //  ifndef __enthSystem_hh
