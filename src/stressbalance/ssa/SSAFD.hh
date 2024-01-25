// Copyright (C) 2004--2024 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <array>

#include "pism/stressbalance/ssa/SSA.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Viewer.hh"
#include "pism/util/petscwrappers/KSP.hh"
#include "pism/util/petscwrappers/Mat.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {
namespace stressbalance {

//! PISM's SSA solver: the finite difference implementation.
class SSAFD : public SSA {
public:
  SSAFD(std::shared_ptr<const Grid> g);
  virtual ~SSAFD() = default;

  const array::Staggered &integrated_viscosity() const;

  const array::Vector &driving_stress() const;

  void compute_residual(const Inputs &inputs, const array::Vector &velocity, array::Vector &result);

protected:

  // Re-implemented by SSAFD_Regional
  virtual void init_impl();

  DiagnosticList diagnostics_impl() const;

  void pc_setup_bjacobi();

  void pc_setup_asm();

  // Re-implemented by SSAFD_Regional
  virtual void solve(const Inputs &inputs);

  void picard_iteration(const Inputs &inputs, double nuH_regularization,
                        double nuH_iter_failure_underrelax);

  void picard_manager(const Inputs &inputs, double nuH_regularization,
                      double nuH_iter_failure_underrelax);

  void picard_strategy_regularization(const Inputs &inputs);

  void compute_average_ice_hardness(const array::Scalar1 &thickness, const array::Array3D &enthalpy,
                                    const array::CellType1 &cell_type, array::Staggered &result);

  void compute_nuH(const array::Scalar1 &ice_thickness, const array::CellType2 &cell_type,
                   const pism::Vector2d *const *velocity, const array::Staggered &hardness,
                   double nuH_regularization, array::Staggered1 &result);

  void compute_nuH_everywhere(const array::Scalar1 &ice_thickness,
                              const pism::Vector2d *const *velocity,
                              const array::Staggered &hardness, double nuH_regularization,
                              array::Staggered &result);

  void compute_nuH_cfbc(const array::Scalar1 &ice_thickness,
                        const array::CellType2 &cell_type,
                        const pism::Vector2d* const* velocity,
                        const array::Staggered &hardness, double nuH_regularization,
                        array::Staggered &result);

  std::array<double, 2> compute_nuH_norm(const array::Staggered &nuH,
                                         array::Staggered &nuH_old);

  // Re-implemented by SSAFD_Regional
  virtual void compute_driving_stress(const array::Scalar &ice_thickness,
                                      const array::Scalar1 &surface_elevation,
                                      const array::CellType1 &cell_type,
                                      const array::Scalar1 *no_model_mask,
                                      array::Vector &result) const;

  void initialize_iterations(const Inputs &inputs);

  void fd_operator(const Inputs &inputs, const pism::Vector2d *const *velocity,
                   const array::Staggered1 &nuH, const array::CellType1 &cell_type, Mat *A,
                   array::Vector *Ax);

  void assemble_matrix(const Inputs &inputs, const array::Vector1 &velocity,
                       const array::Staggered1 &nuH, const array::CellType1 &cell_type, Mat A);

  void assemble_rhs(const Inputs &inputs, const array::CellType1 &cell_type,
                    const array::Vector &driving_stress, array::Vector &result);

  void write_system_petsc(const std::string &namepart);

  void update_nuH_viewers(const array::Staggered &nuH);

  void fracture_induced_softening(const array::Scalar &fracture_density,
                                  array::Staggered &ice_hardness);

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
  // temprary storage used to compute the nuH term
  array::Array2D<Work> m_work;

  array::CellType2 m_cell_type;

  petsc::KSP m_KSP;
  petsc::Mat m_A;
  array::Vector m_rhs;            // right hand side
  array::Vector m_taud;

  array::Vector1 m_velocity_old;
  const double m_scaling;

  // product of the FD matrix and the current guess
  array::Vector m_Ax;

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
SSA * SSAFDFactory(std::shared_ptr<const Grid> grid);

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSAFD_H_ */
