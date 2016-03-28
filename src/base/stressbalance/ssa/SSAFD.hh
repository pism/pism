// Copyright (C) 2004--2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/Viewer.hh"
#include "base/util/petscwrappers/KSP.hh"
#include "base/util/petscwrappers/Mat.hh"

namespace pism {
namespace stressbalance {

//! PISM's SSA solver: the finite difference implementation.
class SSAFD : public SSA
{
  friend class SSAFD_nuH;
public:
  SSAFD(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e);
  virtual ~SSAFD();

  virtual void update(bool fast, double sea_level, const IceModelVec2S &melange_back_pressure);

protected:
  virtual void init_impl();

  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict);

  virtual void pc_setup_bjacobi();

  virtual void pc_setup_asm();
  
  virtual void solve();

  virtual void picard_iteration(double nuH_regularization,
                                double nuH_iter_failure_underrelax);

  virtual void picard_manager(double nuH_regularization,
                              double nuH_iter_failure_underrelax);

  virtual void picard_strategy_regularization();

  virtual void compute_hardav_staggered();

  virtual void compute_nuH_staggered(IceModelVec2Stag &result,
                                     double epsilon);

  virtual void compute_nuH_staggered_cfbc(IceModelVec2Stag &result,
                                          double nuH_regularization);

  virtual void compute_nuH_norm(double &norm,
                                double &norm_change);

  virtual void assemble_matrix(bool include_basal_shear, Mat A);

  virtual void assemble_rhs();

  virtual void write_system_petsc(const std::string &namepart);

  virtual void write_system_matlab(const std::string &namepart);
private:
  PetscErrorCode write_system_matlab_c(const petsc::Viewer &viewer,
                                       const std::string &file_name,
                                       const std::string &cmdstr,
                                       double year);
protected:
  virtual void update_nuH_viewers();

  virtual void set_diagonal_matrix_entry(Mat A, int i, int j,
                                         double value);

  virtual bool is_marginal(int i, int j, bool ssa_dirichlet_bc);

  virtual void fracture_induced_softening();

  // objects used internally
  IceModelVec2Stag hardness, nuH, nuH_old;
  IceModelVec2 m_work;
  petsc::KSP m_KSP;
  petsc::Mat m_A;
  IceModelVec2V m_b;            // right hand side
  double m_scaling;

  const IceModelVec2S *fracture_density, *m_melange_back_pressure;
  IceModelVec2V m_velocity_old;

  unsigned int m_default_pc_failure_count,
    m_default_pc_failure_max_count;
  
  bool view_nuh;
  petsc::Viewer::Ptr nuh_viewer;
  int nuh_viewer_size;

  bool dump_system_matlab;

  class KSPFailure : public RuntimeError {
  public:
    KSPFailure(const char* reason);
  };

  class PicardFailure : public RuntimeError {
  public:
    PicardFailure(const std::string &message);
  };
};

//! Constructs a new SSAFD
SSA * SSAFDFactory(IceGrid::ConstPtr , EnthalpyConverter::Ptr);

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSAFD_H_ */
