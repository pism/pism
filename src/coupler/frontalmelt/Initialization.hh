/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#ifndef POINITIALIZATION_H
#define POINITIALIZATION_H

#include "pism/coupler/FrontalMeltModel.hh"

namespace pism {
namespace frontalmelt {

/*! Frontal melt model "modifier" that helps with initialization.
 *
 * This modifier saves *all* fields a frontal melt model provides as a part of the model state and re-loads
 * them during initialization so that they are available *before* the first time step in a
 * re-started run.
 *
 * It is
 *
 * - not visible to the user,
 * - is added automatically, and
 * - does not have a corresponding "keyword" in frontalmelt::Factory.
 */
class InitializationHelper : public FrontalMeltModel {
public:
  InitializationHelper(IceGrid::ConstPtr g, std::shared_ptr<FrontalMeltModel> in);

private:
  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);

  const IceModelVec2S& frontal_melt_rate_impl() const;

  IceModelVec2S::Ptr m_frontal_melt_rate;
};

} // end of namespace frontalmelt
} // end of namespace pism

#endif /* PFMINITIALIZATION_H */
