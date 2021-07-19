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
#include <mpi.h>

#include "VariableMetadata.hh"

namespace pism {

class File;
class Logger;

//! \brief A general class for reading and accessing time-series.
class Timeseries {
public:
  Timeseries(units::System::Ptr units_system, const std::string &name);

  std::string name() const;
  
  void read(const File &nc, const std::string &time_units, const Logger &log);

  /*!
   * Return the value at time `time`, using linear interpolation
   */
  double operator()(double time) const;
  /*!
   * Return `j`-th value.
   */
  double operator[](unsigned int j) const;

  /*!
   * Return the average over the time interval `[t, t + dt]` estimated using `N`
   * sub-intervals.
   */
  double average(double t, double dt, unsigned int N) const;

  VariableMetadata& variable();

private:
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
