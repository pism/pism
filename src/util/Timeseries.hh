// Copyright (C) 2009, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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

#ifndef __Timeseries_hh
#define __Timeseries_hh

#include <deque>
#include <mpi.h>

#include "VariableMetadata.hh"

namespace pism {

class IceGrid;
class PIO;
class Time;
class Logger;

//! \brief A general class for reading and accessing time-series.
/*!
  \section timeseries_overview Scalar time-series

  This class provides random access to time-series values. It is used to
  implement forcing with scalar time-dependent parameters (such as
  paleo-climate forcing).

  Note that every processor stores the whole time-series and calling append()
  repeatedly will use a lot of memory.

  Please use DiagnosticTimeseries to output long time-series.

  \subsection timeseries_example An example

  The following snippet from PAForcing::init() illustrates creating a Timeseries
  object and reading data from a file.

  \code
  delta_T = new Timeseries(grid.com, grid.rank, "delta_T", "time");
  ierr = delta_T->set_string("units", "Kelvin", ""); CHKERRQ(ierr);
  ierr = delta_T->set_dimension_units("years", ""); CHKERRQ(ierr);
  ierr = delta_T->set_attr("long_name", "near-surface air temperature offsets");
  CHKERRQ(ierr);
  
  ierr = delta_T->read(dT_file); CHKERRQ(ierr);
  \endcode

  Call
  \code
  double offset = (*delta_T)(time);
  \endcode
  to get the value corresponding to the time "time", in this case in years. The
  value returned will be computed using linear interpolation.

  It is also possible to get an n-th value from a time-series: just use square brackets:
  \code
  double offset = (*delta_T)[10];
  \endcode
*/
class Timeseries {
public:
  Timeseries(const IceGrid &g, const std::string &name, const std::string &dimension_name);
  Timeseries(MPI_Comm com, units::System::Ptr units_system,
             const std::string &name, const std::string &dimension_name);
  
  void read(const PIO &nc, const Time &time_manager, const Logger &log);
  void write(const PIO &nc) const;
  double operator()(double time) const;
  double operator[](unsigned int j) const;
  double average(double t, double dt, unsigned int N) const;
  void append(double value, double a, double b);

  void reset();

  TimeseriesMetadata& variable();
  TimeseriesMetadata& dimension();
  TimeBoundsMetadata& bounds();

  const TimeseriesMetadata& variable() const;
  const TimeseriesMetadata& dimension() const;
  const TimeBoundsMetadata& bounds() const;

  const std::vector<double> &times() const;
  const std::vector<double> &time_bounds() const;
  const std::vector<double> &values() const;

  void scale(double scaling_factor);

  std::string name() const;

  bool get_use_bounds() const;
  void set_use_bounds(bool flag);

private:
  MPI_Comm m_com;

  bool m_use_bounds;

  TimeseriesMetadata m_dimension;
  TimeseriesMetadata m_variable;
  TimeBoundsMetadata m_bounds;

  std::vector<double> m_time;
  std::vector<double> m_values;
  std::vector<double> m_time_bounds;

  void set_bounds_units();
  void private_constructor(MPI_Comm com, const std::string &dimension_name);
  void report_range(const Logger &log);
};

} // end of namespace pism

#endif // __Timeseries_hh
