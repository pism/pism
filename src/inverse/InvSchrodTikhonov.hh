// Copyright (C) 2012  David Maxwell
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

#ifndef INVSCHRODTIKHONOV_HH_LD0HAMOE
#define INVSCHRODTIKHONOV_HH_LD0HAMOE

#include "SNESCallbacks.hh"
#include "IceModelVec.hh"
#include "FETools.hh"

class InvSchrodTikhonov {
public:
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  InvSchrodTikhonov( IceGrid  &grid, IceModelVec2V &f);
  ~InvSchrodTikhonov();

  void set_c(IceModelVec2S &c) {
    m_c = &c;
  }

  void setDirichletData(IceModelVec2Int &locs, IceModelVec2V &values) {
    m_dirichletLocations = &locs;
    m_dirichletValues = &values;
  }

  PetscErrorCode solve(bool &success);

  IceModelVec2V &solution() {
    return m_u;
  }

  SNESConvergedReason reason() {
    return m_reason;
  }
  
  std::string reasonDescription() {
    return std::string(SNESConvergedReasons[m_reason]);
  }
  
  PetscErrorCode assembleFunction( DMDALocalInfo *info, PISMVector2 **x, PISMVector2 **f);
  PetscErrorCode assembleJacobian( DMDALocalInfo *info, PISMVector2 **x, Mat J);
  
protected:

  PetscErrorCode construct();
  PetscErrorCode destruct();

  IceGrid &m_grid;
  IceModelVec2S   *m_c;
  IceModelVec2V   *m_f;
  
  IceModelVec2Int *m_dirichletLocations;
  IceModelVec2V   *m_dirichletValues;
  PetscReal        m_dirichletWeight;
  
  IceModelVec2V  m_uGlobal;
  IceModelVec2V  m_u;
  IceModelVec2V  m_r;
  
  SNESDMCallbacks<InvSchrodTikhonov,PISMVector2 **> m_callbacks;

  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;
  
  SNES m_snes;
  DM   m_da;
  Mat  m_J;
  
  SNESConvergedReason m_reason;
/*  
  IceModelVec2V &solution();
  
  PetscErrorCode linearizeAt( IceModelVec2S &c);
  
  void setDirichletData( IceModelVec2Int &locations,
                         IceModelVec2V &velocities );
  
  PetscErrorCode solve();
  
  PetscErrorCode evalObjective(IceModelVec2S &dc, PetscReal *OUTPUT);
  PetscErrorCode evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT);

  PetscErrorCode evalGradObjective(IceModelVec2S &dc, IceModelVec2S &gradient);
  PetscErrorCode evalGradPenaltyReduced(IceModelVec2V &du, IceModelVec2S &gradient);
*/
  // PetscErrorCode evalGradPenalty(IceModelVec2V &du, IceModelVec2V &gradient);

};


#endif /* end of include guard: INVSCHRODTIKHONOV_HH_LD0HAMOE */
