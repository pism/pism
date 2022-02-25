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
#include "pism/util/IceModelVec2S.hh"
#include "pism/util/IceModelVec2Stag.hh"
#include "pism/util/IceModelVec2S.hh"
#include "pism/util/IceModelVec2V.hh"

namespace pism {

class Geometry;

namespace hydrology {

class EmptyingProblem : public Component {
public:
  EmptyingProblem(IceGrid::ConstPtr g);
  virtual ~EmptyingProblem() = default;

  void update(const Geometry &geometry,
              const IceModelVec2S *no_model_mask,
              const IceModelVec2S &water_input_rate,
              bool recompute_potential = true);

  // output
  const IceModelVec2V& flux() const;

  // diagnostics
  const IceModelVec2S& remaining_water_thickness() const;
  const IceModelVec2V& effective_water_velocity() const;
  const IceModelVec2S& potential() const;
  const IceModelVec2S& adjustment() const;
  const IceModelVec2S& sinks() const;

  DiagnosticList diagnostics() const;

protected:

  virtual void compute_raw_potential(const IceModelVec2S &ice_thickness,
                                     const IceModelVec2S &ice_bottom_surface,
                                     IceModelVec2S &result) const;

  void compute_potential(const IceModelVec2S &ice_thickness,
                         const IceModelVec2S &ice_bottom_surface,
                         const IceModelVec2S &domain_mask,
                         IceModelVec2S &result);

  void compute_velocity(const IceModelVec2S &hydraulic_potential,
                        const IceModelVec2S &mask,
                        IceModelVec2Stag &result) const;

  void compute_mask(const array::CellType0 &cell_type,
                    const IceModelVec2S *no_model_mask,
                    IceModelVec2S &result) const;

  Array2SGhosted<1> m_potential;
  IceModelVec2S m_tmp;
  IceModelVec2S m_bottom_surface;
  Array2SGhosted<1> m_W;
  IceModelVec2Stag m_Vstag;
  IceModelVec2Stag m_Qsum;
  Array2SGhosted<1> m_domain_mask;

  IceModelVec2V m_Q;
  IceModelVec2V m_q_sg;
  IceModelVec2S m_adjustment;
  IceModelVec2S m_sinks;

  double m_dx;
  double m_dy;

  double m_eps_gradient;
  double m_speed;
  double m_tau;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* EMPTYINGPROBLEM_H */
