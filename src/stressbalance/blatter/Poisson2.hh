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
#ifndef PISM_POISSON2_H
#define PISM_POISSON2_H

#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/iceModelVec3Custom.hh"
#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {
namespace stressbalance {

class Poisson2 : public ShallowStressBalance {
public:
  Poisson2(IceGrid::ConstPtr grid);
  virtual ~Poisson2();

  void update(const Inputs &inputs, bool);

  const IceModelVec2S& solution() const;
  const IceModelVec2S& exact() const;

  double error() const;
protected:
  void exact_solution(IceModelVec2S &result);

  IceModelVec2S m_solution;
  IceModelVec2S m_exact;

  petsc::DM m_da;
  petsc::SNES m_snes;

  struct CallbackData {
    DM da;
    Poisson2 *solver;
  };

  CallbackData m_callback_data;

  void compute_residual(DMDALocalInfo *info, const double **xg, double **yg);
  static PetscErrorCode function_callback(DMDALocalInfo *info, const double **x, double **f,
                                          CallbackData *data);
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_POISSON2_H */
