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

#ifndef _INVSSAFORWARDPROBLEM_H_
#define _INVSSAFORWARDPROBLEM_H_

#include "SSAFEM.hh"
#include "InvTaucParameterization.hh"

//! \file 
//! \brief Class for implementing the hard parts of a 'siple' 
// NonlinearForwardProblem for the SSA.
/*!\file
Discussion goes here about what the forward problem is, what its
linearization (T) is, what the adjoint of said linearization (T^*)
is, and the inner products.
*/

//! Forward problem for the map from yeild stress to velocities in the SSA
class InvSSAForwardProblem : public SSAFEM
{

public:

  InvSSAForwardProblem(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, 
    EnthalpyConverter &e, InvTaucParameterization &tp,
                 const NCConfigVariable &c) 
           : SSAFEM(g,b,i,e,c), 
             m_KSP(0), m_KSP_B(0), m_MatA(0), m_MatB(0),
             m_VecU(0), m_VecZ2(0),
             m_VecZ(0), m_VecRHS2(0),
             m_VecV(0), m_VecRHS(0),
             m_l2_weight(NULL),
             m_tauc_param(tp),
             m_reassemble_T_matrix_needed(true),
             m_forward_F_needed(true)
  {
    PetscErrorCode ierr = allocate_ksp();
    if (ierr != 0) {
      PetscPrintf(grid.com, "FATAL ERROR: InvSSAForwardProblem allocation failed.\n");
      PISMEnd();
    }    
    ierr = allocate_store();
    if (ierr != 0) {
      PetscPrintf(grid.com, "FATAL ERROR: InvSSAForwardProblem allocation failed.\n");
      PISMEnd();
    }    
  };

  virtual ~InvSSAForwardProblem()
  {
    deallocate_store();
    deallocate_ksp();
  }

  virtual PetscErrorCode init(PISMVars &vars);

  PetscErrorCode set_initial_velocity_guess(IceModelVec2V &v);

  PetscErrorCode set_tauc(IceModelVec2S &tauc );

  PetscErrorCode setup_vars();

  PetscErrorCode solveF(IceModelVec2V &result);

  PetscErrorCode solveT( IceModelVec2S &d, IceModelVec2V &result);

  PetscErrorCode solveTStar( IceModelVec2V &r, IceModelVec2S &result);

  PetscErrorCode domainIP(IceModelVec2S &a, IceModelVec2S &b, PetscScalar *OUTPUT);

  PetscErrorCode rangeIP(IceModelVec2V &a, IceModelVec2V &b, PetscScalar *OUTPUT);

  PetscErrorCode domainIP(Vec a, Vec b, PetscScalar *OUTPUT);

  PetscErrorCode rangeIP(Vec a, Vec b, PetscScalar *OUTPUT);

protected:
  PetscErrorCode allocate_ksp();  

  PetscErrorCode deallocate_ksp();

  PetscErrorCode allocate_store();  

  PetscErrorCode deallocate_store();

  PetscErrorCode domainIP_core(PetscReal **A, PetscReal**B, PetscScalar *OUTPUT);

  PetscErrorCode rangeIP_core(PISMVector2 **A, PISMVector2**B, PetscScalar *OUTPUT);
  
  PetscErrorCode solveF_core();

  PetscErrorCode assemble_T_matrix();

  PetscErrorCode assemble_DomainNorm_matrix();

  PetscErrorCode assemble_T_rhs( PISMVector2 **gvel, PetscReal **gdtau, PISMVector2 **grhs);

  PetscErrorCode assemble_TStarA_rhs( PISMVector2 **R, PISMVector2 **RHS);

  PetscErrorCode assemble_TStarB_rhs( PISMVector2 **Z, PISMVector2 **U, PetscScalar **RHS );

  PetscErrorCode compute_range_l2_area(PetscScalar *OUTPUT);
  
  // PetscErrorCode assemble_TStar_rhs();

  KSP m_KSP, m_KSP_B;
  //! Matrices involved in T and T*
  Mat m_MatA, m_MatB;
  Vec m_VecU;
  //! Left- and right-hand side for linear vector problems.
  Vec m_VecZ2, m_VecZ, m_VecRHS2;
  //! Left- and right-hand side for linear scalar problems.
  Vec m_VecV, m_VecRHS;

  // Optional weight for a weighted L2 norm in the range.
  IceModelVec2S *m_l2_weight;

  // Store for values of dtauc_dxi at the quad points.
  
  PetscReal *m_dtauc_dp_store;

  PetscReal m_range_l2_area;

  InvTaucParameterization &m_tauc_param;

  bool m_reassemble_T_matrix_needed, m_forward_F_needed;
};

//_INVSSAFORWARDPROBLEM_H_
#endif 
