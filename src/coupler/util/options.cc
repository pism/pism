/* Copyright (C) 2018 PISM Authors
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

#include "options.hh"

#include "pism/util/Context.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Component.hh"

namespace pism {

ForcingOptions::ForcingOptions(const Context &ctx,
                               const std::string &prefix) {

  const Logger &log = *ctx.log();
  const Config &config = *ctx.config();

  {
    auto file = config.get_string(prefix + ".file");

    if (not file.empty()) {
      this->filename = file;
      log.message(2,
                  "  - Reading boundary conditions from '%s'...\n",
                  file.c_str());
    } else {
      this->filename = process_input_options(ctx.com(), ctx.config()).filename;

      log.message(2,
                  "  - Option %s_file is not set. Trying the input file '%s'...\n",
                  prefix.c_str(), this->filename.c_str());
    }

    if (this->filename.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "PISM ERROR: %s.file is empty and no input ('-i') file found.",
                                    prefix.c_str());
    }
  }

  {
    auto P = config.get_double(prefix + ".period");

    if (P < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "invalid %s_period %d (period length cannot be negative)",
                                    prefix.c_str(), (int)P);
    }
    this->period = P;
  }

  {
    auto ref_year = config.get_double(prefix + ".reference_year");

    this->reference_time = units::convert(ctx.unit_system(), ref_year, "years", "seconds");
  }
}

} // end of namespace pism
