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

#include "error_handling.hh"
#include <petsc.h>

#include <stdexcept>
#include <cstdarg>


namespace pism {

RuntimeError::RuntimeError(const std::string &message)
  : std::runtime_error(message) {
  // empty
}


RuntimeError RuntimeError::formatted(const char format[], ...) {
  char buffer[8192];
  va_list argp;

  va_start(argp, format);
  vsnprintf(buffer, sizeof(buffer), format, argp);
  va_end(argp);

  return RuntimeError(buffer);
}

RuntimeError::~RuntimeError() throw() {
  // empty
}

void RuntimeError::add_context(const std::string &message) {
  m_context.push_back(message);
}

void RuntimeError::add_context(const char format[], ...) {
  char buffer[8192];
  va_list argp;

  va_start(argp, format);
  vsnprintf(buffer, sizeof(buffer), format, argp);
  va_end(argp);

  // convert to std::string to avoid recursion
  this->add_context(std::string(buffer));
}

void RuntimeError::print(MPI_Comm com) {
  PetscErrorCode ierr = 0;
  std::string error = "PISM ERROR: ",
    message = this->what();

  std::string padding = std::string(error.size(), ' ');

  // replace newlines with newlines plus padding
  size_t k = message.find("\n", 0);
  while (k != std::string::npos) {
    message.insert(k+1, padding);
    k = message.find("\n", k+1);
  }

  // print the error message with "PISM ERROR:" in front:
  ierr = PetscPrintf(com,
                     "%s%s\n", error.c_str(), message.c_str()); CHKERRCONTINUE(ierr);

  // compute how much padding we need to align things:
  std::string while_str = std::string(error.size(), ' ') + "while ";
  padding = std::string(while_str.size() + 1, ' '); // 1 extra space

  // loop over "context" messages
  std::vector<std::string>::const_iterator j = m_context.begin();
  while (j != m_context.end()) {
    message = *j;

    // replace newlines with newlines plus padding
    k = message.find("\n", 0);
    while (k != std::string::npos) {
      message.insert(k+1, padding);
      k = message.find("\n", k+1);
    }

    // print a "context" message
    ierr = PetscPrintf(com,
                       "%s%s\n", while_str.c_str(), message.c_str()); CHKERRCONTINUE(ierr);
    ++j;
  }
}

/** Handle fatal PISM errors by printing an informative error message.
 *
 * (Since these are fatal there is nothing else that can be done.)
 *
 * Should be called from a catch(...) block *only*.
 */
void handle_fatal_errors(MPI_Comm com) {
  PetscErrorCode ierr;
  try {
    throw;                      // re-throw the current exception
  }
  catch (RuntimeError &e) {
    e.print(com);
  }
  catch (std::exception &e) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
                       "\n"
                       "PISM ERROR: Caught a C++ standard library exception: \"%s\".\n"
                       "            This is probably a bug in PISM.\n"
                       "            Please send a report to uaf-pism@alaska.edu\n"
                       "\n",
                       e.what()); CHKERRCONTINUE(ierr);
  } catch (...) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
                       "\n"
                       "PISM ERROR: Caught an unexpected exception.\n"
                       "            This is probably a bug in PISM.\n"
                       "            Please send a report to uaf-pism@alaska.edu\n"
                       "\n");
    CHKERRCONTINUE(ierr);
  }
}

void check_c_call(int errcode, int success,
                  const char* function_name, const char *file, int line) {
  if (errcode != success) {
    throw RuntimeError::formatted("External library function %s failed at %s:%d",
                                  function_name, file, line);
  }
}

void check_petsc_call(int errcode,
                      const char* function_name, const char *file, int line) {
  // tell PETSc to print the error message
  CHKERRCONTINUE(errcode);
  check_c_call(errcode, 0, function_name, file, line);
}

ParallelSection::ParallelSection(MPI_Comm com)
  : m_failed(false), m_com(com) {
  // empty
}

ParallelSection::~ParallelSection() {
  // empty
}

//! @brief Indicates a failure of a parallel section.
/*!
 * This should be called from a `catch (...) { ... }` block **only**.
 */
void ParallelSection::failed() {
  int rank = 0;
  MPI_Comm_rank(m_com, &rank);

  PetscPrintf(MPI_COMM_SELF, "PISM ERROR: ### Rank %d message:\n", rank);

  handle_fatal_errors(MPI_COMM_SELF);

  PetscPrintf(MPI_COMM_SELF, "PISM ERROR: ### Rank %d message ends here.\n", rank);

  m_failed = true;
}

void ParallelSection::reset() {
  m_failed = false;
}

void ParallelSection::check() {
  int success_flag = m_failed ? 0 : 1;
  int success_flag_global = 0;

  MPI_Allreduce(&success_flag, &success_flag_global, 1, MPI_INT, MPI_LAND, m_com);

  if (success_flag_global == 0) {
    throw RuntimeError("Failure in a parallel section. See error messages above for more.");
  }
}

} // end of namespace pism
