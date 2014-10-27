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

#ifndef _ERROR_HANDLING_H_
#define _ERROR_HANDLING_H_

#include <mpi.h>                // MPI_Comm
#include <stdexcept>
#include <string>
#include <vector>

namespace pism {

class RuntimeError : public std::runtime_error {
public:
  RuntimeError(const std::string &message);
  ~RuntimeError() throw();

  //! @brief build a RuntimeError with a formatted message
  static RuntimeError formatted(const char format[], ...);

  //! @brief Add a message providing some context. This way we can (sort of)
  //! get a stack trace even though C++ exceptions do not provide one.
  void add_context(const std::string &message);
  void add_context(const char format[], ...);
  std::vector<std::string> get_context() const;
protected:
  std::vector<std::string> m_context;
};

void handle_fatal_errors(MPI_Comm com);

void check_c_call(int errcode, int success, const char* function_name,
                  const char *file, int line);

void check_petsc_call(int errcode, const char* function_name,
                      const char *file, int line);

#define PISM_CHK(errcode,success,name) do { pism::check_c_call(errcode, success, name, __FILE__, __LINE__); } while (0)
#define PISM_PETSC_CHK(errcode,name) do { pism::check_petsc_call(errcode, name, __FILE__, __LINE__); } while (0)

} // end of namespace pism

#endif /* _ERROR_HANDLING_H_ */
