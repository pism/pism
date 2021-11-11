// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2019, 2021 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>

#include "pism_options.hh"
#include "pism_utilities.hh"
#include "VariableMetadata.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"
#include "Logger.hh"

namespace pism {

//! \brief Print a usage message.
void show_usage(const Logger &log, const std::string &execname, const std::string &usage) {
  log.message(1,
             "%s is a PISM (http://www.pism.io) executable.\n"
             "Options cheat-sheet:\n",
             execname.c_str());
  log.message(1, usage);
  log.message(1,
             "Parallel run using N processes (typical case):  mpiexec -n N %s ...\n"
             "For more help with PISM:\n"
             "  1. download PDF User's Manual or read the online version on the website:\n"
             "       http://www.pism.io\n"
             "  2. read browser for technical details:\n"
             "       http://www.pism.io/doxygen\n"
             "  3. view issues/bugs at source host: https://github.com/pism/pism/issues\n"
             "  4. do '%s -help | grep foo' to see PISM and PETSc options with 'foo'.\n"
             "  5. email for help:  uaf-pism@alaska.edu\n",
             execname.c_str(), execname.c_str());
}

//! @brief In a single call a driver program can provide a usage string to
//! the user and check if required options are given, and if not, end.
bool show_usage_check_req_opts(const Logger &log,
                               const std::string &execname,
                               const std::vector<std::string> &required_options,
                               const std::string &usage) {
  const bool
    keep_running = false,
    terminate = true;

  log.message(2, "%s %s\n", execname.c_str(), pism::revision);

  if (options::Bool("-version", "stop after printing print PISM version")) {
    log.message(2, pism::version());
    return terminate;
  }

  if (options::Bool("-usage", "print PISM usage")) {
    show_usage(log, execname, usage);
    return terminate;
  }

  // go through list of required options, and if not given, fail
  bool req_absent = false;
  for (auto opt : required_options) {
    if (not options::Bool(opt, "a required option")) {
      req_absent = true;
      log.error("PISM ERROR: option %s required\n", opt.c_str());
    }
  }

  if (req_absent) {
    log.error("\n");
    show_usage(log, execname, usage);
    return terminate;
  }

  // show usage message with -help, but don't stop
  if (options::Bool("-help", "print help on all options")) {
    show_usage(log, execname, usage);
  }
  return keep_running;
}

} // end of namespace pism
