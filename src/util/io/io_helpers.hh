/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#ifndef _IO_HELPERS_H_
#define _IO_HELPERS_H_

#include <string>
#include <vector>
#include <mpi.h>

#include "IO_Flags.hh"
#include "pism/util/Units.hh"
#include "pism/util/interpolation.hh"

namespace pism {

class VariableMetadata;
class SpatialVariableMetadata;
class IceGrid;
class File;
class Time;
class Logger;
class Context;
class Config;

namespace io {

void regrid_spatial_variable(SpatialVariableMetadata &var,
                             const IceGrid& grid, const File &nc,
                             RegriddingFlag flag, bool do_report_range,
                             bool allow_extrapolation,
                             double default_value,
                             InterpolationType type,
                             double *output);

void regrid_spatial_variable(SpatialVariableMetadata &var,
                             const IceGrid& grid, const File &nc,
                             unsigned int t_start,
                             RegriddingFlag flag, bool do_report_range,
                             bool allow_extrapolation,
                             double default_value,
                             InterpolationType type,
                             double *output);

void read_spatial_variable(const SpatialVariableMetadata &var,
                           const IceGrid& grid, const File &nc,
                           unsigned int time, double *output);

void write_spatial_variable(const SpatialVariableMetadata &var,
                            const IceGrid& grid, const File &nc,
                            const double *input);

void define_dimension(const File &nc, unsigned long int length,
                      const VariableMetadata &metadata);

void define_time(const File &file, const Context &ctx);

void define_time(const File &nc, const std::string &name, const std::string &calendar,
                 const std::string &units, units::System::Ptr unit_system);

void append_time(const File &file, const Config &ctx, double time_seconds);
void append_time(const File &nc, const std::string &name, double time_seconds);

void define_spatial_variable(const SpatialVariableMetadata &var,
                             const IceGrid &grid, const File &nc,
                             IO_Type nctype);

void define_timeseries(const VariableMetadata& var,
                       const std::string &dimension_name,
                       const File &nc, IO_Type nctype);

void define_time_bounds(const VariableMetadata& var,
                        const std::string &dimension_name,
                        const std::string &bounds_name,
                        const File &nc, IO_Type nctype = PISM_DOUBLE);

void read_timeseries(const File &nc, const VariableMetadata &metadata,
                     const Logger &log, std::vector<double> &data);

void write_timeseries(const File &nc, const VariableMetadata &metadata,
                      size_t t_start, const std::vector<double> &data);

void read_time_bounds(const File &nc,
                      const VariableMetadata &metadata,
                      const Logger &log, std::vector<double> &data);

void write_time_bounds(const File &nc, const VariableMetadata &metadata,
                       size_t t_start, const std::vector<double> &data);

void read_time_info(const Logger &log,
                    std::shared_ptr<units::System> unit_system,
                    const File &file,
                    const std::string &time_name,
                    const std::string &time_units,
                    std::vector<double> &times,
                    std::vector<double> &bounds);

std::string time_dimension(units::System::Ptr unit_system,
                           const File &file,
                           const std::string &variable_name);

void read_attributes(const File &nc, const std::string &variable_name, VariableMetadata &variable);

void write_attributes(const File &nc, const VariableMetadata &variable, IO_Type nctype);

void read_valid_range(const File &nc, const std::string &name, VariableMetadata &variable);

bool file_exists(MPI_Comm com, const std::string &filename);

void move_if_exists(MPI_Comm com, const std::string &file_to_move, int rank_to_use = 0);

void remove_if_exists(MPI_Comm com, const std::string &file_to_remove, int rank_to_use = 0);

} // end of namespace io
} // end of namespace pism

#endif /* _IO_HELPERS_H_ */
