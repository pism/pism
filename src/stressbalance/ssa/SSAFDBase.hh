/* Copyright (C) 2024 PISM Authors
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

#ifndef PISM_SSAFDBASE_H
#define PISM_SSAFDBASE_H

#include <petscmat.h>

#include "pism/stressbalance/ssa/SSA.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {
namespace stressbalance {

/*!
 * A base class containing the FD discretization of the SSA system.
 *
 * Does *not* include any implementation details related to non-linear iterations.
 */
class SSAFDBase : public SSA {
public:
  SSAFDBase(std::shared_ptr<const Grid> g, bool regional_mode);

  const array::Staggered &integrated_viscosity() const;

  const array::Vector &driving_stress() const;

  void compute_residual(const Inputs &inputs, const array::Vector &velocity, array::Vector &result);
protected:
  void initialize_iterations(const Inputs &inputs);

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

  void compute_driving_stress(const array::Scalar &ice_thickness,
                              const array::Scalar1 &surface_elevation,
                              const array::CellType1 &cell_type,
                              const array::Scalar1 *no_model_mask, const EnthalpyConverter &EC,
                              array::Vector &result) const;

  void adjust_driving_stress(const array::Scalar &ice_thickness,
                             const array::Scalar1 &surface_elevation,
                             const array::CellType1 &cell_type, const array::Scalar1 *no_model_mask,
                             array::Vector &driving_stress) const;

  void compute_average_ice_hardness(const array::Scalar1 &thickness, const array::Array3D &enthalpy,
                                    const array::CellType1 &cell_type, array::Staggered &result) const;

  void assemble_rhs(const Inputs &inputs, const array::CellType1 &cell_type,
                    const array::Vector &driving_stress, double bc_scaling, array::Vector &result) const;

  void fd_operator(const Geometry &geometry, const array::Scalar *bc_mask, double bc_scaling,
                   const array::Scalar &basal_yield_stress,
                   IceBasalResistancePlasticLaw *basal_sliding_law,
                   const pism::Vector2d *const *velocity, const array::Staggered1 &nuH,
                   const array::CellType1 &cell_type, Mat *A, array::Vector *Ax) const;

  void fracture_induced_softening(const array::Scalar1 &fracture_density,
                                  double n_glen,
                                  array::Staggered &ice_hardness);

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

  //! temprary storage used to compute the nuH term (ghosted, but ghost values are
  //! computed "redundantly" and not communicated)
  array::Array2D<Work> m_work;

  //! ice hardness
  array::Staggered m_hardness;

  //! viscosity times thickness
  array::Staggered1 m_nuH;

  array::CellType2 m_cell_type;

  //! right hand side
  array::Vector m_rhs;

  //! driving stress
  array::Vector m_taud;

  //! scaling used for diagonal matrix elements at Dirichlet BC locations
  const double m_bc_scaling;

  //! if true, the driving stress is adjusted to avoid issues at domain boundaries in
  //! "regional" model configurations
  const bool m_regional_mode;
};

}
} // namespace pism

#endif /* PISM_SSAFDBASE_H */
