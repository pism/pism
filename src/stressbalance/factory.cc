/* Copyright (C) 2017, 2018, 2020, 2021, 2023, 2024 PISM Authors
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

#include <memory>

#include "pism/regional/SIAFD_Regional.hh"
#include "pism/stressbalance/SSB_Modifier.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/WeertmanSliding.hh"
#include "pism/stressbalance/blatter/Blatter.hh"
#include "pism/stressbalance/blatter/BlatterMod.hh"
#include "pism/stressbalance/sia/SIAFD.hh"
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFD_SNES.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "pism/util/Context.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace stressbalance {

std::shared_ptr<StressBalance> create(const std::string &model,
                                      std::shared_ptr<const Grid> grid,
                                      bool regional) {

  auto config = grid->ctx()->config();

  if (model == "blatter") {
    int Mz = config->get_number("stress_balance.blatter.Mz");
    int C = config->get_number("stress_balance.blatter.coarsening_factor");

    auto blatter = std::make_shared<Blatter>(grid, Mz, C);
    auto mod = std::make_shared<BlatterMod>(blatter);

    return std::make_shared<StressBalance>(grid, blatter, mod);
  }

  auto ssa_method = config->get_string("stress_balance.ssa.method");
  std::shared_ptr<ShallowStressBalance> sliding;
  if (member(model, {"none", "sia"})) {
    sliding = std::make_shared<ZeroSliding>(grid);
  } else if (member(model, {"prescribed_sliding", "prescribed_sliding+sia"})) {
    sliding = std::make_shared<PrescribedSliding>(grid);
  } else if (member(model, {"weertman_sliding", "weertman_sliding+sia"})) {
    sliding = std::make_shared<WeertmanSliding>(grid);
  } else if (member(model, {"ssa", "ssa+sia"})) {
    if (ssa_method == "fd") {
      sliding = std::make_shared<SSAFD>(grid, regional);
    } else if (ssa_method == "fd_snes") {
      sliding = std::make_shared<SSAFD_SNES>(grid, regional);
    } else {
      sliding = std::make_shared<SSAFEM>(grid);
    }
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid stress balance model: %s", model.c_str());
  }

  std::shared_ptr<SSB_Modifier> modifier;
  if (member(model, {"none", "ssa", "prescribed_sliding", "weertman_sliding"})) {
    modifier = std::make_shared<ConstantInColumn>(grid);
  } else if (member(model, {"prescribed_sliding+sia", "weertman_sliding+sia", "ssa+sia", "sia"})) {
    if (regional) {
      modifier = std::make_shared<SIAFD_Regional>(grid);
    } else {
      modifier = std::make_shared<SIAFD>(grid);
    }
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid stress balance model: %s", model.c_str());
  }

  return std::make_shared<StressBalance>(grid, sliding, modifier);
}

} // end of namespace stressbalance
} // end of namespace pism
