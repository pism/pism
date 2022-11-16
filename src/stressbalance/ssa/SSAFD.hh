// Copyright (C) 2004--2022 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Viewer.hh"
#include "pism/util/petscwrappers/KSP.hh"
#include "pism/util/petscwrappers/Mat.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {
namespace stressbalance {

//! PISM's SSA solver: the finite difference implementation.
class SSAFD : public SSA
{
public:
  SSAFD(IceGrid::ConstPtr g);
  virtual ~SSAFD() = default;

  const array::Staggered & integrated_viscosity() const;
protected:
  virtual void init_impl();

  virtual DiagnosticList diagnostics_impl() const;

  virtual void pc_setup_bjacobi();

  virtual void pc_setup_asm();
  
  virtual void solve(const Inputs &inputs);

  virtual void picard_iteration(const Inputs &inputs,
                                double nuH_regularization,
                                double nuH_iter_failure_underrelax);

  virtual void picard_manager(const Inputs &inputs,
                              double nuH_regularization,
                              double nuH_iter_failure_underrelax);

  virtual void picard_strategy_regularization(const Inputs &inputs);

  virtual void compute_hardav_staggered(const Inputs &inputs);

  virtual void compute_nuH_staggered(const Geometry &geometry,
                                     double nuH_regularization,
                                     array::Staggered &result);

  virtual void compute_nuH_staggered_cfbc(const Geometry &geometry,
                                          double nuH_regularization,
                                          array::Staggered &result);

  virtual void compute_nuH_norm(double &norm,
                                double &norm_change);

  virtual void assemble_matrix(const Inputs &inputs,
                               bool include_basal_shear, Mat A);

  virtual void assemble_rhs(const Inputs &inputs);

  virtual void write_system_petsc(const std::string &namepart);

  virtual void update_nuH_viewers();

  virtual bool is_marginal(int i, int j, bool ssa_dirichlet_bc);

  virtual void fracture_induced_softening(const array::Scalar *fracture_density);

  // objects used internally
  array::Staggered m_hardness;
  array::Staggered1 m_nuH, m_nuH_old;

  struct Work {
    // u_x on the i offset
    double u_x;
    // v_x on the i offset
    double v_x;
    // weight for the i offset
    double w_i;
    // u_y on the j offset
    double u_y;
    // v_y on the j offset
    double v_y;
    // weight for the j offset
    double w_j;
  };
  array::Array2D<Work> m_work;

  petsc::KSP m_KSP;
  petsc::Mat m_A;
  array::Vector m_b;            // right hand side

  array::Vector1 m_velocity_old;
  const double m_scaling;

  unsigned int m_default_pc_failure_count,
    m_default_pc_failure_max_count;
  
  bool m_view_nuh;
  std::shared_ptr<petsc::Viewer> m_nuh_viewer;
  int m_nuh_viewer_size;

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
SSA * SSAFDFactory(IceGrid::ConstPtr grid);

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSAFD_H_ */
