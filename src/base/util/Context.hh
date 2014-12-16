/* Copyright (C) 2014 PISM Authors
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

#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include <string>
#include <mpi.h>
#include "PISMConfig.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "PISMUnits.hh"

namespace pism {

class Context {
public:
  Context(Config &config);
  MPI_Comm com() const;

  Config &config();
  const Config& config() const;

  Vars& variables();
  const Vars& variables() const;

  double convert(double value, const std::string &unit1, const std::string &unit2) const;
  UnitSystem unit_system() const;

  Time::Ptr time();
  const Time::Ptr time() const;

  Profiling profiling() const;
private:
  Config m_config;
  Vars m_variables;
  mutable Profiling m_profiling;
  Time::Ptr m_time;

  // Hide copy constructor / assignment operator.
  Context(Context const &);
  Context & operator=(Context const &);
};

} // end of namespace pism

#endif /* _CONTEXT_H_ */
