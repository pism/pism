/* Copyright (C) 2014, 2015 PISM Authors
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

#include "pism_memory.hh"

namespace pism {

namespace units {
class System;
}

class Config;
class EnthalpyConverter;
class Time;
class Profiling;
class Logger;

class Context {
public:
  typedef PISM_SHARED_PTR(units::System) UnitsSystemPtr;
  typedef PISM_SHARED_PTR(Config) ConfigPtr;
  typedef PISM_SHARED_PTR(const Config) ConstConfigPtr;
  typedef PISM_SHARED_PTR(EnthalpyConverter) EnthalpyConverterPtr;
  typedef PISM_SHARED_PTR(Time) TimePtr;
  typedef PISM_SHARED_PTR(const Time) ConstTimePtr;
  typedef PISM_SHARED_PTR(Logger) LoggerPtr;
  typedef PISM_SHARED_PTR(const Logger) ConstLoggerPtr;

  typedef PISM_SHARED_PTR(Context) Ptr;
  typedef PISM_SHARED_PTR(const Context) ConstPtr;

  Context(MPI_Comm c, UnitsSystemPtr sys,
          ConfigPtr conf, EnthalpyConverterPtr EC, TimePtr t,
          LoggerPtr log,
          const std::string &p);

  MPI_Comm com() const;
  UnitsSystemPtr unit_system() const;
  ConstConfigPtr config() const;
  EnthalpyConverterPtr enthalpy_converter() const;
  ConstTimePtr time() const;
  const std::string& prefix() const;
  const Profiling& profiling() const;

  ConstLoggerPtr log() const;
  LoggerPtr log();

  ConfigPtr config();
  TimePtr time();
private:
  class Impl;
  PISM_SHARED_PTR(Impl) m_impl;
  // disable copying and assignments
  Context(const Context& other);
  Context & operator=(const Context &);
};

//! Create a default context using options.
Context::Ptr context_from_options(MPI_Comm com, const std::string &prefix);

} // end of namespace pism

#endif /* _CONTEXT_H_ */
