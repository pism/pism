/* Copyright (C) 2015, 2016, 2017, 2018 PISM Authors
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

namespace pism {

class VariableMetadata;
class SpatialVariableMetadata;
class TimeseriesMetadata;
class TimeBoundsMetadata;
class IceGrid;
class PIO;
class Time;
class Logger;
class Context;
class Config;

namespace io {

void regrid_spatial_variable(SpatialVariableMetadata &var,
                             const IceGrid& grid, const PIO &nc,
                             RegriddingFlag flag, bool do_report_range,
                             bool allow_extrapolation,
                             double default_value,
                             InterpolationType type,
                             double *output);

void regrid_spatial_variable(SpatialVariableMetadata &var,
                             const IceGrid& grid, const PIO &nc,
                             unsigned int t_start,
                             RegriddingFlag flag, bool do_report_range,
                             bool allow_extrapolation,
                             double default_value,
                             InterpolationType type,
                             double *output);

void read_spatial_variable(const SpatialVariableMetadata &var,
                           const IceGrid& grid, const PIO &nc,
                           unsigned int time, double *output);

void write_spatial_variable(const SpatialVariableMetadata &var,
                            const IceGrid& grid, const PIO &nc,
                            const double *input);

void define_dimension(const PIO &nc, unsigned long int length,
                      const VariableMetadata &metadata);

void define_time(const PIO &file, const Context &ctx);

void define_time(const PIO &nc, const std::string &name, const std::string &calendar,
                 const std::string &units, units::System::Ptr unit_system);

void append_time(const PIO &file, const Config &ctx, double time_seconds);
void append_time(const PIO &nc, const std::string &name, double time_seconds);

void define_spatial_variable(const SpatialVariableMetadata &var,
                             const IceGrid &grid, const PIO &nc,
                             IO_Type nctype,
                             const std::string &variable_order);

void define_timeseries(const TimeseriesMetadata& var,
                       const PIO &nc, IO_Type nctype);

void define_time_bounds(const TimeBoundsMetadata& var,
                        const PIO &nc, IO_Type nctype);

void read_timeseries(const PIO &nc, const TimeseriesMetadata &metadata,
                     const Time &time, const Logger &log, std::vector<double> &data);

void write_timeseries(const PIO &nc, const TimeseriesMetadata &metadata,
                      size_t t_start, const std::vector<double> &data,
                      IO_Type nctype = PISM_DOUBLE);

void write_timeseries(const PIO &nc, const TimeseriesMetadata &metadata,
                      size_t t_start, double data,
                      IO_Type nctype = PISM_DOUBLE);

void read_time_bounds(const PIO &nc,
                      const TimeBoundsMetadata &metadata,
                      const Time &time, const Logger &log, std::vector<double> &data);

void write_time_bounds(const PIO &nc, const TimeBoundsMetadata &metadata,
                       size_t t_start, const std::vector<double> &data,
                       IO_Type nctype = PISM_DOUBLE);

void read_attributes(const PIO &nc, const std::string &variable_name, VariableMetadata &variable);

void write_attributes(const PIO &nc, const VariableMetadata &variable, IO_Type nctype);

void read_valid_range(const PIO &nc, const std::string &name, VariableMetadata &variable);

bool file_exists(MPI_Comm com, const std::string &filename);

} // end of namespace io
} // end of namespace pism

#endif /* _IO_HELPERS_H_ */
