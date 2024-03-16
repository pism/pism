// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023, 2024 PISM Authors
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

#include "ClimateIndexWeights.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"

namespace pism {

ClimateIndexWeights::ClimateIndexWeights(const Context &ctx) {
  m_index.reset(
      new ScalarForcing(ctx, "climate_index", "climate_index", "1", "1", "climate index"));

  const Config &config = *ctx.config();

  m_ref_value = config.get_number("climate_index.ref_value");
  m_max_value = config.get_number("climate_index.max_value");
}

std::array<double,3> ClimateIndexWeights::update_weights(double t, double dt) {

  double CI = m_index->value(t + 0.5 * dt);

  // compute interpolation weights in update
  double W0 = 1.0 - std::min(CI, m_ref_value) / m_ref_value;

  double W1 = (std::max(CI, m_ref_value) - m_ref_value) / (1.0 - m_ref_value);
  W1 = std::min(W1, 1.0);

  double W1X = (std::max(CI, 1.0) - 1.0) / (m_max_value - 1.0);

  return {W0, W1, W1X};
}

} // end of namespace pism
