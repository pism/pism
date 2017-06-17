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

#include <unistd.h>
#include <cstdarg>

#include "Logger.hh"
#include "pism_const.hh"
#include "pism_options.hh"

namespace pism {

struct Logger::Impl {
  MPI_Comm com;
  bool enabled;
  int threshold;
};

Logger::Logger(MPI_Comm com, int threshold)
  : m_impl(new Impl) {

  m_impl->com = com;
  m_impl->enabled = true;
  m_impl->threshold = threshold;
}

Logger::~Logger() {
  delete m_impl;
}

void Logger::message(int threshold, const char format[], ...) const {
  if ((not m_impl->enabled) or threshold > m_impl->threshold) {
    return;
  }

  char buffer[8192];
  va_list argp;

  va_start(argp, format);
  vsnprintf(buffer, sizeof(buffer), format, argp);
  va_end(argp);

  message_impl(buffer);
}

void Logger::message(int threshold, const std::string &buffer) const {
  if ((not m_impl->enabled) or threshold > m_impl->threshold) {
    return;
  }

  message_impl(buffer.c_str());
}

void Logger::message_impl(const char buffer[]) const {
  verbPrintf(1, m_impl->com, buffer);
}

void Logger::set_threshold(int level) {
  m_impl->threshold = level;
}
int Logger::get_threshold() const {
  return m_impl->threshold;
}
void Logger::disable() {
  m_impl->enabled = false;
}

void Logger::enable() {
  m_impl->enabled = true;
}

Logger::Ptr logger_from_options(MPI_Comm com) {
  Logger::Ptr result(new Logger(com, 2));

  options::Integer verbosity("-verbose", "set logger verbosity threshold",
                             getVerbosityLevel());

  setVerbosityLevel(verbosity);

  result->set_threshold(verbosity);

  return result;
}

} // end of namespace pism
