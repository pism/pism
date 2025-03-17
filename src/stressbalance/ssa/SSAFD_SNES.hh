/* Copyright (C) 2024, 2025 PISM Authors
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

#ifndef PISM_SSAFD_SNES_H
#define PISM_SSAFD_SNES_H

#include "pism/stressbalance/ssa/SSAFDBase.hh"

#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/petscwrappers/DM.hh"
// #include "pism/util/petscwrappers/Vec.hh"

namespace pism {
namespace stressbalance {

class SSAFD_SNES : public SSAFDBase {
public:
  SSAFD_SNES(std::shared_ptr<const Grid> grid, bool regional_mode);

  void solve(const Inputs &inputs);

  double tolerance() const;

  const array::Vector &residual() const;
private:
  DiagnosticList diagnostics_impl() const;

  //! residual (diagnostic)
  array::Vector m_residual;

  petsc::SNES m_snes;
  std::shared_ptr<petsc::DM> m_DA;

  struct CallbackData {
    DM da;
    SSAFD_SNES *solver;
    const Inputs *inputs;
  };

  CallbackData m_callback_data;

  void compute_jacobian(const Inputs &inputs, Vector2d const *const * velocity, Mat J);

  static PetscErrorCode function_callback(DMDALocalInfo *info,
                                          Vector2d const *const * velocity,
                                          Vector2d **result,
                                          CallbackData *);
  static PetscErrorCode jacobian_callback(DMDALocalInfo *info,
                                          Vector2d const *const * velocity,
                                          Mat A, Mat J,
                                          CallbackData *data);
};

} // namespace stressbalance
} // namespace pism

#endif /* PISM_SSAFD_SNES_H */
