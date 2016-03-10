// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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
#include <cassert>

#include "pism_options.hh"
#include "VariableMetadata.hh"
#include "PISMConfigInterface.hh"

#include "error_handling.hh"
#include "Logger.hh"

namespace pism {

//! Determine verbosity level from user options.
/*!
\verbatim
   level  option        meaning
   -----  ------        -------
   0      -verbose 0    never print to std out AT ALL!
   1      -verbose 1    less verbose than default: thresh must be 1 to print
   2      -verbose 2    DEFAULT
   3      -verbose 3    somewhat verbose
          -verbose      same as "-verbose 3"
   4      -verbose 4    fairly verbose
   5      -verbose 5    very verbose: print everything
\endverbatim
See verbPrintf().
 */
void verbosityLevelFromOptions() {

  setVerbosityLevel(2);

  options::String just_verbose("-verbose", "verbosity level", "");

  if (just_verbose.is_set() and just_verbose->empty()) {
    setVerbosityLevel(3);
  } else {
    options::Integer verbose("-verbose", "verbosity level", 2);
    if (verbose.is_set()) {
      setVerbosityLevel(verbose);
    }
  }
}

//! \brief Print a usage message.
void show_usage(const Logger &log, const std::string &execname, const std::string &usage) {
  log.message(1,
             "%s is a PISM (http://www.pism-docs.org) executable.\n"
             "Options cheat-sheet:\n",
             execname.c_str());
  log.message(1, usage);
  log.message(1,
             "Parallel run using N processes (typical case):  mpiexec -n N %s ...\n"
             "For more help with PISM:\n"
             "  1. download PDF User's Manual:\n"
             "       http://www.pism-docs.org/wiki/lib/exe/fetch.php?media=pism_manual.pdf\n"
             "  2. read browser for technical details:\n"
             "       http://www.pism-docs.org/doxy/html/index.html\n"
             "  3. view issues/bugs at source host: https://github.com/pism/pism/issues\n"
             "  4. do '%s -help | grep foo' to see PISM and PETSc options with 'foo'.\n"
             "  5. email for help:  help@pism-docs.org\n",
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

  if (options::Bool("-usage", "print PISM usage")) {
    show_usage(log, execname, usage);
    return terminate;
  }

  // go through list of required options, and if not given, fail
  bool req_absent = false;
  for (size_t k = 0; k < required_options.size(); ++k) {
    if (not options::Bool(required_options[k], "a required option")) {
      req_absent = true;
      log.error("PISM ERROR: option %s required\n", required_options[k].c_str());
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
