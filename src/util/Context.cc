/* Copyright (C) 2014, 2015, 2017, 2019, 2021, 2023, 2024, 2025 PISM Authors
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

#include "pism/util/Context.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/Units.hh"
#include "pism/util/Time.hh"
#include "pism/util/Logger.hh"
#include "pism/util/EnthalpyConverter.hh"
#include <memory>

namespace pism {

class Context::Impl {
public:
  Impl(MPI_Comm c,
       std::shared_ptr<units::System> sys,
       std::shared_ptr<Config> conf,
       std::shared_ptr<EnthalpyConverter> EC,
       std::shared_ptr<Time> t,
       std::shared_ptr<Logger> log,
       const std::string &p)
    : com(c), unit_system(sys), config(conf), enthalpy_converter(EC), time(t), prefix(p),
      logger(log) {
    // empty
  }
  MPI_Comm com;
  std::shared_ptr<units::System> unit_system;
  std::shared_ptr<Config> config;
  std::shared_ptr<EnthalpyConverter> enthalpy_converter;
  std::shared_ptr<Time> time;
  std::string prefix;
  Profiling profiling;
  std::shared_ptr<Logger> logger;
};

Context::Context(MPI_Comm c, std::shared_ptr<units::System> sys,
                 std::shared_ptr<Config> config, std::shared_ptr<EnthalpyConverter> EC, std::shared_ptr<Time> t,
                 std::shared_ptr<Logger> L,
                 const std::string &p)
  : m_impl(new Impl(c, sys, config, EC, t, L, p)) {
  // empty
}

Context::~Context() {
  delete m_impl;
}

MPI_Comm Context::com() const {
  return m_impl->com;
}

int Context::size() const {
  int S = 0;
  MPI_Comm_size(m_impl->com, &S);
  return S;
}

int Context::rank() const {
  int R = 0;
  MPI_Comm_rank(m_impl->com, &R);
  return R;
}

std::shared_ptr<units::System> Context::unit_system() const {
  return m_impl->unit_system;
}

std::shared_ptr<Config> Context::config() {
  return m_impl->config;
}

std::shared_ptr<const Config> Context::config() const {
  return m_impl->config;
}

std::shared_ptr<EnthalpyConverter> Context::enthalpy_converter() const {
  return m_impl->enthalpy_converter;
}

std::shared_ptr<Time> Context::time() {
  return m_impl->time;
}

std::shared_ptr<const Time> Context::time() const {
  return m_impl->time;
}

const std::string& Context::prefix() const {
  return m_impl->prefix;
}

const Profiling& Context::profiling() const {
  return m_impl->profiling;
}

std::shared_ptr<const Logger> Context::log() const {
  return m_impl->logger;
}

std::shared_ptr<Logger> Context::log() {
  return m_impl->logger;
}

std::shared_ptr<Context> context_from_options(MPI_Comm com,
                                              const std::string &prefix,
                                              bool print) {
  // unit system
  auto sys = std::make_shared<units::System>();

  // logger
  auto logger = std::make_shared<Logger>(com, 1);

  // configuration parameters
  auto config = config_from_options(com, sys);

  logger->set_threshold(static_cast<int>(config->get_number("output.runtime.verbosity")));

  if (print) {
    print_config(*logger, 3, *config);
  }

  // time manager
  auto time = std::make_shared<Time>(com, config, *logger, sys);

  // enthalpy converter
  auto EC = std::make_shared<EnthalpyConverter>(*config);

  return std::make_shared<Context>(com, sys, config, EC, time, logger, prefix);
}


} // end of namespace pism
