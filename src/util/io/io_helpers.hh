/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023, 2024, 2025 PISM Authors
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

#ifndef PISM_IO_HELPERS_H
#define PISM_IO_HELPERS_H

#include "OutputWriter.hh"
#include <memory>
#include <string>
#include <vector>
#include <mpi.h>

namespace pism {

namespace units {
class System;
}

class VariableMetadata;
class Grid;
class File;
class Logger;

namespace grid {
class InputGridInfo;
class DistributedGridInfo;
}

class LocalInterpCtx;

class OutputFile;

namespace io {

std::string time_dimension(std::shared_ptr<units::System> unit_system,
                           const File &file,
                           const std::string &variable_name);

void check_input_grid(const grid::InputGridInfo &input_grid,
                      const grid::DistributedGridInfo& internal_grid,
                      const std::vector<double> &internal_z_levels,
                      bool allow_extrapolation);

void regrid_spatial_variable(const VariableMetadata &variable,
                             const Grid& internal_grid,
                             const LocalInterpCtx &lic,
                             const File &file,
                             const Logger &log,
                             double *output);

void read_spatial_variable(const VariableMetadata &variable,
                           const Grid& grid, const File &file,
                           unsigned int time, double *output);

std::vector<double> read_1d_variable(const File &file, const std::string &name,
                                     const std::string &units,
                                     std::shared_ptr<units::System> unit_system);

std::vector<double> read_timeseries_variable(const File &file, const std::string &variable_name,
                                             const std::string &units,
                                             std::shared_ptr<units::System> unit_system,
                                             size_t start, size_t count);

std::vector<double> read_bounds(const File &file, const std::string &bounds_variable_name,
                                const std::string &units,
                                std::shared_ptr<units::System> unit_system);

void read_time_info(std::shared_ptr<units::System> unit_system, const File &file,
                    const std::string &time_name, const std::string &time_units,
                    std::vector<double> &times, std::vector<double> &bounds);

VariableMetadata read_attributes(const File &file, const std::string &variable_name,
                                 std::shared_ptr<units::System> unit_system);

// writing utilities

/*!
 * Define the time dimension
 */
void define_time(const OutputFile &output_file, const VariableMetadata &metadata,
                 bool with_bounds);

void move_if_exists(MPI_Comm com, const std::string &file_to_move, int rank_to_use = 0);

void remove_if_exists(MPI_Comm com, const std::string &file_to_remove, int rank_to_use = 0);

} // end of namespace io
} // end of namespace pism

#endif /* PISM_IO_HELPERS_H */
