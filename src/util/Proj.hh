/* Copyright (C) 2016, 2017, 2019, 2020, 2024 PISM Authors
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

#ifndef PROJ_WRAPPER_H
#define PROJ_WRAPPER_H

#include <string>

#include <proj.h>

#include "pism/util/error_handling.hh"

namespace pism {

//! A wrapper for PJ that makes sure `pj_destroy` is called.
class Proj {
public:
  Proj(const std::string &input, const std::string &output) {
    // add +type=crs if it is not there yet
    std::string input_with_crs(input);
    if ((input.find("+proj") != std::string::npos or input.find("+init") != std::string::npos) and
        input.find("+type=crs") == std::string::npos) {
      input_with_crs += " +type=crs";
    }

    m_pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX, input_with_crs.c_str(), output.c_str(), 0);

    if (m_pj == 0) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "Failed to initialize projection transformation '%s' to '%s' (errno: %d, %s).",
          input_with_crs.c_str(), output.c_str(), proj_errno(0), proj_errno_string(proj_errno(0)));
    }
  }

  ~Proj() {
    proj_destroy(m_pj);
  }
  operator PJ*() {
    return m_pj;
  }
private:
  PJ *m_pj;
};

} // end of namespace pism


#endif /* PROJ_WRAPPER_H */
