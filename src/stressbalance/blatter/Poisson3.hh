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
#ifndef PISM_POISSON3_H
#define PISM_POISSON3_H

#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/iceModelVec3Custom.hh"
#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {
namespace stressbalance {

class Poisson3 : public ShallowStressBalance {
public:
  Poisson3(IceGrid::ConstPtr grid);
  virtual ~Poisson3();

  void update(const Inputs &inputs, bool);

  IceModelVec3Custom::Ptr solution() const;
  IceModelVec3Custom::Ptr exact() const;
protected:
  void exact_solution(double b, double H, IceModelVec3Custom &result);

  IceModelVec3Custom::Ptr m_solution;
  IceModelVec3Custom::Ptr m_exact;

  petsc::DM m_da;
  petsc::Vec m_x;
  petsc::SNES m_snes;

  struct CallbackData {
    DM da;
    Poisson3 *solver;
  };

  CallbackData m_callback_data;

  void compute_residual(DMDALocalInfo *info, const double ***xg, double ***yg);
  static PetscErrorCode function_callback(DMDALocalInfo *info, const double ***x, double ***f,
                                          CallbackData *data);
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_POISSON3_H */
