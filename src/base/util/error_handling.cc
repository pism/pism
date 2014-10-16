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

#include "error_handling.hh"
#include <petsc.h>

#include <stdexcept>

namespace pism {

RuntimeError::RuntimeError(const std::string &message)
  : std::runtime_error(message) {
  // empty
}

RuntimeError::~RuntimeError() throw() {
  // empty
}

void RuntimeError::add_context(const std::string &message) {
  m_context.push_back(message);
}

std::vector<std::string> RuntimeError::get_context() const {
  return m_context;
}

/** Handle fatal PISM errors by printing an informative error message.
 *
 * (Since these are fatal there is nothing else that can be done.)
 *
 * Should be called from a catch(...) block *only*.
 */
void handle_fatal_errors(MPI_Comm com) {
  try {
    throw;                      // re-throw the current exception
  }
  catch (RuntimeError &e) {
    std::vector<std::string> context = e.get_context();
    std::string error = "PISM ERROR: ",
      message = e.what();

    std::string padding = std::string(error.size(), ' ');

    // replace newlines with newlines plus padding
    size_t k = message.find("\n", 0);
    while (k != std::string::npos) {
      message.insert(k+1, padding);
      k = message.find("\n", k+1);
    }

    // print the error message with "PISM ERROR:" in front:
    PetscPrintf(com,
                "%s%s\n", error.c_str(), message.c_str());

    // compute how much padding we need to align things:
    std::string while_str = std::string(error.size(), ' ') + "while ";
    padding = std::string(while_str.size() + 1, ' '); // 1 extra space

    // loop over "context" messages
    std::vector<std::string>::const_iterator j = context.begin();
    while (j != context.end()) {
      message = *j;

      // replace newlines with newlines plus padding
      k = message.find("\n", 0);
      while (k != std::string::npos) {
        message.insert(k+1, padding);
        k = message.find("\n", k+1);
      }

      // print a "context" message
      PetscPrintf(com,
                  "%s%s\n", while_str.c_str(), message.c_str());
      ++j;
    }
  }
  catch (std::exception &e) {
    PetscPrintf(PETSC_COMM_SELF,
                "PISM ERROR: caught a C++ standard library exception: %s.\n"
                "     This is probably a bug in PISM."
                " Please send a report to help@pism-docs.org\n",
                e.what());
  }
  catch (...) {
    PetscPrintf(PETSC_COMM_SELF,
                "PISM ERROR: caught an unexpected exception.\n"
                "     This is probably a bug in PISM."
                " Please send a report to help@pism-docs.org\n");
  }
}

} // end of namespace pism
