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

#ifndef INVSSATIKHONOV_HH_NSO6ITYB
#define INVSSATIKHONOV_HH_NSO6ITYB

#include "SSAFEM.hh"
#include "InvTaucParameterization.hh"
#include "Functional.hh"
#include <memory> // for std::auto_ptr


class InvSSATikhonov : public SSAFEM
{

public:

  InvSSATikhonov(IceGrid &g, IceBasalResistancePlasticLaw &b,
    EnthalpyConverter &e, InvTaucParameterization &tp,
    const NCConfigVariable &c);

  virtual ~InvSSATikhonov();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode set_tauc_fixed_locations(IceModelVec2Int &locations)
  { 
    m_fixed_tauc_locations = &locations;
    return 0;
  }

  PetscErrorCode set_zeta( IceModelVec2S &zeta);

  PetscErrorCode linearizeAt( IceModelVec2S &c, bool &success);

  PetscErrorCode evalObjective(IceModelVec2S &dc, PetscReal *OUTPUT);
  PetscErrorCode evalGradObjective(IceModelVec2S &dc, IceModelVec2S &gradient);

  PetscErrorCode evalPenalty(IceModelVec2V &du, PetscReal *OUTPUT);
  PetscErrorCode evalGradPenaltyReduced(IceModelVec2V &dr, IceModelVec2S &gradient);

protected:

  PetscErrorCode construct();
  PetscErrorCode destruct();

  IceGrid &m_grid;
  IceModelVec2S   *m_zeta;

  IceModelVec2Int *m_fixed_tauc_locations;
  IceModelVec2S   *m_misfit_weight;

  InvTaucParameterization &m_tauc_param;
  
  IceModelVec2V  m_vGlobal;
  IceModelVec2V  m_v;
  IceModelVec2V  m_adjointRHS;
  
  FEElementMap m_element_index;
  FEQuadrature m_quadrature;
  FEDOFMap     m_dofmap;

  std::auto_ptr< Functional<IceModelVec2S> > m_designFunctional;
  std::auto_ptr< Functional<IceModelVec2V> > m_penaltyFunctional;
  
  KSP  m_ksp;
  Mat  m_Jadjoint;
  
  SNESConvergedReason m_reason;
};

#endif /* end of include guard: INVSSATIKHONOV_HH_NSO6ITYB */
