/* Copyright (C) 2016, 2017, 2019 PISM Authors
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

#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H
#include <proj_api.h>

#include "pism/util/error_handling.hh"

namespace pism {

//! A wrapper for projPJ that makes sure `pj_free` is called.
class Proj {
public:
  Proj(const std::string &proj_string) {
    m_proj = pj_init_plus(proj_string.c_str());
    if (m_proj == NULL) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Failed to initialize projection '%s'.",
                                    proj_string.c_str());
    }
  }
  ~Proj() {
    pj_free(m_proj);
  }
  operator projPJ() {
    return m_proj;
  }
private:
  projPJ m_proj;
};

} // end of namespace pism


#endif /* PROJ_WRAPPER_H */
