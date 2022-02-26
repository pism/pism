/* Copyright (C) 2019, 2020, 2021, 2022 PISM Authors
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

#ifndef EMPTYINGPROBLEM_H
#define EMPTYINGPROBLEM_H

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/IceModelVec2Stag.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/IceModelVec2V.hh"

namespace pism {

class Geometry;

namespace hydrology {

class EmptyingProblem : public Component {
public:
  EmptyingProblem(IceGrid::ConstPtr g);
  virtual ~EmptyingProblem() = default;

  void update(const Geometry &geometry,
              const array::Scalar *no_model_mask,
              const array::Scalar &water_input_rate,
              bool recompute_potential = true);

  // output
  const IceModelVec2V& flux() const;

  // diagnostics
  const array::Scalar& remaining_water_thickness() const;
  const IceModelVec2V& effective_water_velocity() const;
  const array::Scalar& potential() const;
  const array::Scalar& adjustment() const;
  const array::Scalar& sinks() const;

  DiagnosticList diagnostics() const;

protected:

  virtual void compute_raw_potential(const array::Scalar &ice_thickness,
                                     const array::Scalar &ice_bottom_surface,
                                     array::Scalar &result) const;

  void compute_potential(const array::Scalar &ice_thickness,
                         const array::Scalar &ice_bottom_surface,
                         const array::Scalar &domain_mask,
                         array::Scalar &result);

  void compute_velocity(const array::Scalar &hydraulic_potential,
                        const array::Scalar &mask,
                        IceModelVec2Stag &result) const;

  void compute_mask(const array::CellType0 &cell_type,
                    const array::Scalar *no_model_mask,
                    array::Scalar &result) const;

  array::Scalar1 m_potential;
  array::Scalar m_tmp;
  array::Scalar m_bottom_surface;
  array::Scalar1 m_W;
  IceModelVec2Stag m_Vstag;
  IceModelVec2Stag m_Qsum;
  array::Scalar1 m_domain_mask;

  IceModelVec2V m_Q;
  IceModelVec2V m_q_sg;
  array::Scalar m_adjustment;
  array::Scalar m_sinks;

  double m_dx;
  double m_dy;

  double m_eps_gradient;
  double m_speed;
  double m_tau;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* EMPTYINGPROBLEM_H */
