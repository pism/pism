// Copyright (C) 2009, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2021 Constantine Khroulev
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
#include <mpi.h>                // MPI_Comm

#include "VariableMetadata.hh"

namespace pism {

class IceGrid;
class File;
class Time;
class Logger;

//! \brief A general class for reading and accessing time-series.
/*!
  \section timeseries_overview Scalar time-series

  This class provides random access to time-series values. It is used to
  implement forcing with scalar time-dependent parameters (such as
  paleo-climate forcing).

  \subsection timeseries_example An example

  The following snippet from PAForcing::init() illustrates creating a Timeseries
  object and reading data from a file.

  \code
  Timeseries delta_T(grid, "delta_T", "time");
  delta_T.variable().set_string("units", "Kelvin");

  delta_T->read(dT_file);
  \endcode

  Call
  \code
  double offset = delta_T(time);
  \endcode
  to get the value corresponding to the time "time", in this case in years. The
  value returned will be computed using linear interpolation.

  It is also possible to get an n-th value from a time-series: just use square brackets:
  \code
  double offset = delta_T[10];
  \endcode
*/
class Timeseries {
public:
  Timeseries(const IceGrid &g, const std::string &name);
  Timeseries(MPI_Comm com, units::System::Ptr units_system,
             const std::string &name);

  std::string name() const;
  
  void read(const File &nc, const Time &time_manager, const Logger &log);

  double operator()(double time) const;
  double operator[](unsigned int j) const;
  double average(double t, double dt, unsigned int N) const;

  VariableMetadata& variable();

private:
  MPI_Comm m_com;

  bool m_use_bounds;

  const units::System::Ptr m_unit_system;

  VariableMetadata m_variable;

  std::vector<double> m_time;
  std::vector<double> m_values;
  std::vector<double> m_time_bounds;

  void report_range(const Logger &log);
};

} // end of namespace pism

#endif // __Timeseries_hh
