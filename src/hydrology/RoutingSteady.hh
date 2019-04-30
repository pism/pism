/* Copyright (C) 2019 PISM Authors
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

#ifndef ROUTING_STEADY_STATE_H
#define ROUTING_STEADY_STATE_H

#include "Routing.hh"

namespace pism {
namespace hydrology {

/*!
 * A modification of the "routing" hydrology model that tries to route *all* of the input
 * to the margin (i.e. computing the steady state water distribution).
 */
class RoutingSteady : public Routing {
public:
  RoutingSteady(IceGrid::ConstPtr g);
  virtual ~RoutingSteady();

protected:
  virtual void restart_impl(const PIO &input_file, int record);

  virtual void bootstrap_impl(const PIO &input_file,
                              const IceModelVec2S &ice_thickness);

  virtual void init_impl(const IceModelVec2S &W_till,
                         const IceModelVec2S &W,
                         const IceModelVec2S &P);

  virtual void update_impl(double t, double dt, const Inputs& inputs);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  void compute_velocity(const IceModelVec2Stag &W,
                        const IceModelVec2S &pressure,
                        const IceModelVec2S &bed,
                        const IceModelVec2Int *no_model_mask,
                        IceModelVec2Stag &result) const;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* ROUTING_STEADY_STATE_H */
