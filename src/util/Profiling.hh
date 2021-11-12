/* Copyright (C) 2015 PISM Authors
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

#ifndef PISM_PROFILING_HH
#define PISM_PROFILING_HH

#include <map>
#include <string>
#include <petsclog.h>

namespace pism {

class Profiling {
public:
  Profiling();
  void start() const;
  void report(const std::string &filename) const;
  void begin(const char *name) const;
  void end(const char *name) const;
  void stage_begin(const char *name) const;
  void stage_end(const char *name) const;
private:
  PetscClassId m_classid;
  mutable std::map<std::string, PetscLogEvent> m_events;
  mutable std::map<std::string, PetscLogStage> m_stages;
};

} // end of namespace pism

#endif /* PISM_PROFILING_HH */
