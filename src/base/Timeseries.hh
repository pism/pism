// Copyright (C) 2009 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "NCVariable.hh"
#include "grid.hh"
#include "nc_util.hh"
#include <deque>

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond


//! A general class for reading and accessing time-series.
/*!
  Provides random access to time-series values.

  Note that every processor stores the whole time-series and calling append()
  repeatedly will use a lot of memory.

  Please use DiagnosticTimeseries to output long time-series.
 */
class Timeseries {
public:
  Timeseries(IceGrid * g, string name, string dimension_name);
  Timeseries(MPI_Comm com, PetscMPIInt rank, string name, string dimension_name);
  
  PetscErrorCode read(const char filename[]);
  PetscErrorCode write(const char filename[]);
  double operator()(double time);
  double operator[](unsigned int j) const;
  PetscErrorCode append(double time, double value);
  int length();
  PetscErrorCode set_attr(string name, double value);
  PetscErrorCode set_attr(string name, string value);
  PetscErrorCode set_units(string units, string glaciological_units);
  PetscErrorCode set_dimension_units(string units, string glaciological_units);

protected:
  MPI_Comm com;
  PetscMPIInt rank;
  NCTimeseries dimension, var;

  vector<double> time;
  vector<double> values;
};

//! A class for storing and writing diagnostic time-series.
/*! This version of Timeseries only holds \c buffer_size entries in memory and
  writes to a file every time this limit is exceeded.

  Here is a usage example:

  First, create the DiagnosticTimeseries object and set metadata. This will
  prepare the offsets object to write delta_T(t) time-series to
  pism-delta_T.nc, converting from degrees Celsius (internal units) to degrees
  Kelvin ("glaciological" units). Time will be written in years (%i.e. there is
  no unit conversion there).

  \code
  offsets = new DiagnosticTimeseries(g, "delta_T", "t");
  offsets->set_units("Kelvin", "Celsius");
  offsets->set_dimension_units("years", "");
  offsets->buffer_size = 100; // only store 100 entries; default is 10000
  offsets->set_output_prefix("pism-");
  \endcode

  Once this is set up, one can add calls like

  \code
  offsets->append(t_years, TsOffset);
  \endcode

  to the code. This will store the (t_years, TsOffset) pair, but not
  necessarily cause any I/O.

  Note that every time you exceed the \c buffer_size limit, all the entries are
  written to a file <b> and removed from memory</b>.
 */
class DiagnosticTimeseries : public Timeseries {
public:
  DiagnosticTimeseries(IceGrid * g, string name, string dimension_name);
  DiagnosticTimeseries(MPI_Comm com, PetscMPIInt rank, string name, string dimension_name);
  ~DiagnosticTimeseries();

  PetscErrorCode append(double time, double value);
  PetscErrorCode set_output_prefix(string prefix);
  PetscErrorCode flush();

  size_t buffer_size;
  string output_filename;

protected:
  size_t start;
  deque<double> t, v;
};

#endif // __Timeseries_hh
