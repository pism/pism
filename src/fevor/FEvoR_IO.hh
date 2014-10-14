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

#ifndef _FEVOR_IO_H_
#define _FEVOR_IO_H_

#include <string>
#include <mpi.h>

namespace FEvoR {
class Distribution;
}

namespace pism {
class UnitSystem;
}

int fevor_prepare_file(const std::string &filename,
                       MPI_Comm comm, const pism::UnitSystem &sys,
                       unsigned int n_distributions,
                       unsigned int n_crystals);

int fevor_save_distribution(const std::string &filename,
                            MPI_Comm comm, const pism::UnitSystem &sys,
                            unsigned int index, const FEvoR::Distribution &distribution);

int fevor_load_distribution(const std::string &filename,
                            MPI_Comm comm, const pism::UnitSystem &sys,
                            unsigned int index, FEvoR::Distribution &distribution);

#endif /* _FEVOR_IO_H_ */
