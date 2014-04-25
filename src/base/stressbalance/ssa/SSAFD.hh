// Copyright (C) 2004--2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

namespace pism {

//! PISM's SSA solver: the finite difference implementation.
class SSAFD : public SSA
{
  friend class SSAFD_nuH;
public:
  SSAFD(IceGrid &g, EnthalpyConverter &e, const Config &c);
  virtual ~SSAFD();

  virtual PetscErrorCode init(Vars &vars);

  virtual PetscErrorCode update(bool fast, IceModelVec2S &melange_back_pressure);

  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);
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

  virtual PetscErrorCode compute_hardav_staggered();

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               double epsilon);

  virtual PetscErrorCode compute_nuH_staggered_cfbc(IceModelVec2Stag &result,
                                                    double nuH_regularization);

  virtual PetscErrorCode compute_nuH_norm(double &norm,
                                          double &norm_change);

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs);

  virtual PetscErrorCode write_system_petsc();

  virtual PetscErrorCode write_system_matlab();

  virtual PetscErrorCode update_nuH_viewers();

  virtual PetscErrorCode set_diagonal_matrix_entry(Mat A, int i, int j,
                                                   double value);

  virtual bool is_marginal(int i, int j, bool ssa_dirichlet_bc);

  virtual PetscErrorCode fracture_induced_softening();

  // objects used internally
  IceModelVec2Stag hardness, nuH, nuH_old;
  IceModelVec2 m_work;
  KSP m_KSP;
  Mat m_A;
  Vec m_b;
  double m_scaling;

  IceModelVec2S *fracture_density, *m_melange_back_pressure;

  unsigned int m_default_pc_failure_count,
    m_default_pc_failure_max_count;
  
  bool view_nuh;
  PetscViewer nuh_viewer;
  int nuh_viewer_size;

  bool dump_system_matlab;
};

//! Constructs a new SSAFD
SSA * SSAFDFactory(IceGrid &, EnthalpyConverter &, const Config &);

//! \brief Reports the nuH (viscosity times thickness) product on the staggered
//! grid.
class SSAFD_nuH : public Diag<SSAFD>
{
public:
  SSAFD_nuH(SSAFD *m, IceGrid &g, Vars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

} // end of namespace pism

#endif /* _SSAFD_H_ */
