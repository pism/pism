/* Copyright (C) 2020 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef PISM_BLATTER_H
#define PISM_BLATTER_H

#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/fem/FEM.hh"

#include "grid_hierarchy.hh"    // GridInfo

namespace pism {

namespace fem {
class Element3;
class Q1Element3Face;
} // end of namespace fem

namespace stressbalance {

class Blatter : public ShallowStressBalance {
public:
  Blatter(IceGrid::ConstPtr grid, int Mz, int n_levels, int coarsening_factor);
  virtual ~Blatter();

  void update(const Inputs &inputs, bool);

  IceModelVec3::Ptr velocity_u_sigma() const;
  IceModelVec3::Ptr velocity_v_sigma() const;

protected:
  // u and v components of ice velocity on the sigma grid
  IceModelVec3::Ptr m_u_sigma, m_v_sigma;

  // 3D dof=2 DM used by SNES
  petsc::DM m_da;
  // storage for the solution
  petsc::Vec m_x;

  petsc::SNES m_snes;

  struct CallbackData {
    DM da;
    Blatter *solver;
  };

  CallbackData m_callback_data;
  GridInfo m_grid_info;
  double m_rho_ice_g;
  double m_rho_ocean_g;

  static const int m_Nq = 100;
  static const int m_n_work = 4;

  double m_work[m_n_work][m_Nq];

  Vector2 m_work2[m_n_work][m_Nq];

  void init_impl();

  void define_model_state_impl(const File &output) const;

  void write_model_state_impl(const File &output) const;

  void compute_jacobian(DMDALocalInfo *info, const Vector2 ***x, Mat A, Mat J);

  void jacobian_f(const fem::Element3 &element,
                  const Vector2 *u_nodal,
                  const double *B_nodal,
                  double K[2 * fem::q13d::n_chi][2 * fem::q13d::n_chi]);

  void jacobian_basal(const fem::Q1Element3Face &face,
                      const double *tauc_nodal,
                      const double *f_nodal,
                      const Vector2 *u_nodal,
                      double K[2 * fem::q13d::n_chi][2 * fem::q13d::n_chi]);

  void compute_residual(DMDALocalInfo *info, const Vector2 ***xg, Vector2 ***yg);

  void residual_f(const fem::Element3 &element,
                  const Vector2 *u_nodal,
                  const double *B_nodal,
                  Vector2 *residual);

  void residual_source_term(const fem::Element3 &element,
                            const double *surface,
                            Vector2 *residual);

  void residual_basal(const fem::Element3 &element,
                      const fem::Q1Element3Face &face,
                      const double *tauc_nodal,
                      const double *f_nodal,
                      const Vector2 *u_nodal,
                      Vector2 *residual);

  void residual_lateral(const fem::Element3 &element,
                        const fem::Q1Element3Face &face,
                        const double *z_nodal,
                        const double *sl_nodal,
                        Vector2 *residual);

  static PetscErrorCode jacobian_callback(DMDALocalInfo *info,
                                          const Vector2 ***x,
                                          Mat A, Mat J, CallbackData *data);

  static PetscErrorCode function_callback(DMDALocalInfo *info, const Vector2 ***x, Vector2 ***f,
                                          CallbackData *data);

  void init_2d_parameters(const Inputs &inputs);
  void init_ice_hardness(const Inputs &inputs);

  // Guts of the constructor. This method wraps PETSc calls to simplify error checking.
  PetscErrorCode setup(DM pism_da, int Mz, int n_levels, int coarsening_factor);

  PetscErrorCode setup_2d_storage(DM dm, int dof);

  void set_initial_guess(const IceModelVec3 &u_sigma, const IceModelVec3 &v_sigma);

  void copy_solution();

  void compute_averaged_velocity(IceModelVec2V &result);

  void get_basal_velocity(IceModelVec2V &result);
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_BLATTER_H */
