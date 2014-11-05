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
#include <vector>
#include <mpi.h>

namespace FEvoR {
class Distribution;
}

namespace pism {
class PIO;
class UnitSystem;
}

int fevor_prepare_file(const pism::PIO &nc, const pism::UnitSystem &sys,
                       unsigned int n_distributions,
                       const std::vector<unsigned int> packing_dimensions);

int fevor_save_distribution(const pism::PIO &nc,
                            unsigned int index,
                            unsigned int time_index,
                            const FEvoR::Distribution &distribution);

int fevor_load_distribution(const pism::PIO &nc,
                            unsigned int index,
                            unsigned int time_index,
                            FEvoR::Distribution &distribution);

int fevor_load_particle_positions(const pism::PIO &nc,
                                  unsigned int time_index,
                                  std::vector<double> &x,
                                  std::vector<double> &y,
                                  std::vector<double> &z);

int fevor_save_particle_positions(const pism::PIO &nc,
                                  unsigned int time_index,
                                  std::vector<double> &x,
                                  std::vector<double> &y,
                                  std::vector<double> &z);

int fevor_save_recrystallization_numbers(const pism::PIO &nc,
                                         unsigned int time_index,
                                         std::vector<double> &migration,
                                         std::vector<double> &polygonization);

#endif /* _FEVOR_IO_H_ */
