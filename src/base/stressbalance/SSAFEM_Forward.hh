// Copyright (C) 2009--2011 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#ifndef _SSAFEM_FORWARD_H_
#define _SSAFEM_FORWARD_H_

#include "SSAFEM.hh"

class SSAFEM_Forward : public SSAFEM
{

public:

  SSAFEM_Forward(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, 
                 EnthalpyConverter &e,
                 const NCConfigVariable &c) 
           : SSAFEM(g,b,i,e,c), 
             m_KSP(0), m_MatA(0), m_MatB(0), 
             m_VecZ(0), m_VecRHS2(0),
             m_VecV(0), m_VecRHS(0),
             m_reassemble_T_matrix_needed(true),
             m_forward_F_needed(true)
  {
    allocate_ksp();
  };

  virtual ~SSAFEM_Forward()
  {
    deallocate_ksp();
  }

  PetscErrorCode allocate_ksp();  
  PetscErrorCode deallocate_ksp();

  PetscErrorCode set_initial_velocity_guess(IceModelVec2V &v);

  PetscErrorCode set_tauc(IceModelVec2S &tauc );

  PetscErrorCode setup_vars();

  PetscErrorCode solveF(IceModelVec2V &result);

  PetscErrorCode solveT( IceModelVec2S &d, IceModelVec2V &result);

  PetscErrorCode solveTStar( IceModelVec2V &r, IceModelVec2S &result);

  PetscErrorCode domainIP(IceModelVec2S &a, IceModelVec2S &b, PetscScalar *OUTPUT);

  PetscErrorCode rangeIP(IceModelVec2V &a, IceModelVec2V &b, PetscScalar *OUTPUT);

protected:
  
  PetscErrorCode solveF_core();

  PetscErrorCode assemble_T_matrix();

  PetscErrorCode assemble_DomainNorm_matrix();

  PetscErrorCode assemble_T_rhs( PISMVector2 **gvel, PetscReal **gdtau, PISMVector2 **grhs);

  PetscErrorCode assemble_TStarA_rhs( PISMVector2 **R, PISMVector2 **RHS);

  PetscErrorCode assemble_TStarB_rhs( PISMVector2 **Z, PISMVector2 **U, PetscScalar **RHS );
  
  // PetscErrorCode assemble_TStar_rhs();

  KSP m_KSP, m_KSP_B;
  //! Matrices involved in T and T*
  Mat m_MatA, m_MatB;
  //! Left- and right-hand side for linear vector problems.
  Vec m_VecZ, m_VecRHS2;
  //! Left- and right-hand side for linear scalar problems.
  Vec m_VecV, m_VecRHS;

  bool m_reassemble_T_matrix_needed, m_forward_F_needed;
};

#endif
