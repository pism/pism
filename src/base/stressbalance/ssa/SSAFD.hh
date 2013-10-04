// Copyright (C) 2004--2013 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SSAFD_H_
#define _SSAFD_H_

#include "SSA.hh"
#include <petscksp.h>


//! PISM's SSA solver: the finite difference implementation
class SSAFD : public SSA
{
  friend class SSAFD_nuH;
public:
  SSAFD(IceGrid &g, IceBasalResistancePlasticLaw &b, EnthalpyConverter &e,
        const NCConfigVariable &c) :
    SSA(g,b,e,c)
  {
    PetscErrorCode ierr = allocate_fd();
    if (ierr != 0) {
      PetscPrintf(grid.com, "FATAL ERROR: SSAFD allocation failed.\n");
      PISMEnd();
    }
  }

  virtual ~SSAFD()
  {
    PetscErrorCode ierr = deallocate_fd();
    if (ierr != 0) {
      PetscPrintf(grid.com, "FATAL ERROR: SSAFD de-allocation failed.\n");
      PISMEnd();
    }
  }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict,
                               map<string, PISMTSDiagnostic*> &ts_dict);
protected:
  virtual PetscErrorCode allocate_fd();

  virtual PetscErrorCode deallocate_fd();

  virtual PetscErrorCode pc_setup_bjacobi();

  virtual PetscErrorCode pc_setup_asm();
  
  virtual PetscErrorCode solve();

  virtual PetscErrorCode picard_iteration(unsigned int max_iterations,
                                          double ssa_relative_tolerance,
                                          double nuH_regularization);

  virtual PetscErrorCode strategy_1_regularization();

  virtual PetscErrorCode strategy_2_asm();

  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result);

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               PetscReal epsilon);

  virtual PetscErrorCode compute_nuH_norm(PetscReal &norm,
                                          PetscReal &norm_change);

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs);

  virtual PetscErrorCode write_system_petsc();

  virtual PetscErrorCode write_system_matlab();

  virtual PetscErrorCode update_nuH_viewers();

  virtual PetscErrorCode set_diagonal_matrix_entry(Mat A, int i, int j,
                                                   PetscScalar value);

  virtual bool is_marginal(int i, int j, bool ssa_dirichlet_bc);

  // objects used internally
  IceModelVec2Stag hardness, nuH, nuH_old;
  KSP m_KSP;
  Mat m_A;
  Vec m_b;
  PetscScalar m_scaling;

  unsigned int m_default_pc_failure_count,
    m_default_pc_failure_max_count;
  
  bool view_nuh;
  PetscViewer nuh_viewer;
  PetscInt nuh_viewer_size;

  bool dump_system_matlab;
};

//! Constructs a new SSAFD
SSA * SSAFDFactory(IceGrid &, IceBasalResistancePlasticLaw &,
                   EnthalpyConverter &, const NCConfigVariable &);

//! \brief Reports the nuH (viscosity times thickness) product on the staggered
//! grid.
class SSAFD_nuH : public PISMDiag<SSAFD>
{
public:
  SSAFD_nuH(SSAFD *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _SSAFD_H_ */
