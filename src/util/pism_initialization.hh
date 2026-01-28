/* Copyright (C) 2026 PISM Authors
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

#include <string>
#include <vector>

namespace pism {

/*!
 * Initialize PISM's dependencies (version for use in Python).
 */
void initialize(const char *help);

/*!
 * Initialize PISM's dependencies (version for use in C++).
 */
void initialize(int argc, char **argv, const char *help);

/*!
 * Insert command line options into PETSc's default options database.
 */
void initialize_options(const std::vector<std::string> &args);

/*!
 * Finalize PISM's dependencies.
 */
void finalize();

} // namespace pism
