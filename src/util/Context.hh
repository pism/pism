/* Copyright (C) 2014, 2015, 2016, 2019, 2020, 2021, 2023, 2024 PISM Authors
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

#ifndef PISM_CONTEXT_H
#define PISM_CONTEXT_H

#include <memory>
#include <string>

#include <mpi.h>

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
  Context(MPI_Comm c, std::shared_ptr<units::System> sys, std::shared_ptr<Config> conf,
          std::shared_ptr<EnthalpyConverter> EC, std::shared_ptr<Time> t,
          std::shared_ptr<Logger> log, const std::string &p);
  ~Context();

  MPI_Comm com() const;
  int size() const;
  int rank() const;
  std::shared_ptr<units::System> unit_system() const;
  std::shared_ptr<const Config> config() const;
  std::shared_ptr<EnthalpyConverter> enthalpy_converter() const;
  std::shared_ptr<const Time> time() const;
  const std::string& prefix() const;
  const Profiling& profiling() const;

  std::shared_ptr<const Logger> log() const;
  std::shared_ptr<Logger> log();

  std::shared_ptr<Config> config();
  std::shared_ptr<Time> time();

private:
  class Impl;
  Impl *m_impl;
  // disable copying and assignments
  Context(const Context& other);
  Context & operator=(const Context &);
};

//! Create a default context using options.
std::shared_ptr<Context> context_from_options(MPI_Comm com,
                                              const std::string &prefix,
                                              bool print = false);

} // end of namespace pism

#endif /* PISM_CONTEXT_H */
