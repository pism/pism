/* Copyright (C) 2014, 2015, 2017, 2019, 2021, 2023 PISM Authors
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
#include "pism/util/error_handling.hh"
#include "pism/pism_config.hh"

#if (Pism_USE_PIO==1)
// Why do I need this???
#define _NETCDF
#include <pio.h>
#endif

namespace pism {

class Context::Impl {
public:
  Impl(MPI_Comm c,
       UnitsSystemPtr sys,
       ConfigPtr conf,
       EnthalpyConverterPtr EC,
       TimePtr t,
       LoggerPtr log,
       const std::string &p)
    : com(c), unit_system(sys), config(conf), enthalpy_converter(EC), time(t), prefix(p),
      logger(log), pio_iosys_id(-1) {
    // empty
  }
  MPI_Comm com;
  UnitsSystemPtr unit_system;
  ConfigPtr config;
  EnthalpyConverterPtr enthalpy_converter;
  TimePtr time;
  std::string prefix;
  Profiling profiling;
  LoggerPtr logger;
  int pio_iosys_id;
};

Context::Context(MPI_Comm c, UnitsSystemPtr sys,
                 ConfigPtr config, EnthalpyConverterPtr EC, TimePtr t,
                 LoggerPtr L,
                 const std::string &p)
  : m_impl(new Impl(c, sys, config, EC, t, L, p)) {
  // empty
}

Context::~Context() {

#if (Pism_USE_PIO==1)
  if (m_impl->pio_iosys_id != -1 and
      PIOc_free_iosystem(m_impl->pio_iosys_id) != PIO_NOERR) {
    m_impl->logger->message(1, "Error: failed to de-allocate a ParallelIO I/O system\n");
  }
#endif

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

Context::UnitsSystemPtr Context::unit_system() const {
  return m_impl->unit_system;
}

Context::ConfigPtr Context::config() {
  return m_impl->config;
}

Context::ConstConfigPtr Context::config() const {
  return m_impl->config;
}

Context::EnthalpyConverterPtr Context::enthalpy_converter() const {
  return m_impl->enthalpy_converter;
}

Context::TimePtr Context::time() {
  return m_impl->time;
}

Context::ConstTimePtr Context::time() const {
  return m_impl->time;
}

const std::string& Context::prefix() const {
  return m_impl->prefix;
}

const Profiling& Context::profiling() const {
  return m_impl->profiling;
}

Context::ConstLoggerPtr Context::log() const {
  return m_impl->logger;
}

Context::LoggerPtr Context::log() {
  return m_impl->logger;
}

/*!
 * I/O system id (the ParallelIO library)
 */
int Context::pio_iosys_id() const {
#if (Pism_USE_PIO==1)
  if (m_impl->pio_iosys_id == -1) {
    int ierr = PIOc_set_iosystem_error_handling(PIO_DEFAULT, PIO_BCAST_ERROR, NULL);
    if (ierr != 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Failed to initialize ParallelIO");
    }

    int
      base      = config()->get_number("output.pio.base"),
      stride    = config()->get_number("output.pio.stride"),
      n_writers = config()->get_number("output.pio.n_writers");

    if (n_writers > this->size()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "number of ParallelIO writers (%d)"
                                    " exceeds the communicator size (%d)",
                                    n_writers, this->size());
    }

    ierr = PIOc_Init_Intracomm(m_impl->com, n_writers, stride, base, PIO_REARR_BOX,
                               &m_impl->pio_iosys_id);
    if (ierr != 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Failed to initialize ParallelIO");
    }
  }
#endif
  return m_impl->pio_iosys_id;
}

std::shared_ptr<Context> context_from_options(MPI_Comm com,
                                              const std::string &prefix,
                                              bool print) {
  // unit system
  units::System::Ptr sys(new units::System);

  // logger
  Logger::Ptr logger = logger_from_options(com);

  // configuration parameters
  Config::Ptr config = config_from_options(com, *logger, sys);

  if (print) {
    print_config(*logger, 3, *config);
  }

  // time manager
  Time::Ptr time = std::make_shared<Time>(com, config, *logger, sys);

  // enthalpy converter
  EnthalpyConverter::Ptr EC(new EnthalpyConverter(*config));

  return std::shared_ptr<Context>(new Context(com, sys, config, EC, time, logger, prefix));
}


} // end of namespace pism
