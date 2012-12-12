// Copyright (C) 2009--2012 Constantine Khroulev
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

#include "Timeseries.hh"
#include <algorithm>
#include "pism_const.hh"
#include "IceGrid.hh"
#include "PIO.hh"

Timeseries::Timeseries(IceGrid *g, string name, string dimension_name)
{
  private_constructor(g->com, g->rank, name, dimension_name);
}

Timeseries::Timeseries(MPI_Comm c, PetscMPIInt r, string name, string dimension_name) {
  private_constructor(c, r, name, dimension_name);
}

void Timeseries::private_constructor(MPI_Comm c, PetscMPIInt r, string name, string dimension_name) {
  com = c;
  rank = r;
  dimension.init(dimension_name, dimension_name, com, rank);
  var.init(name, dimension_name, com, rank);
  bounds.init(dimension_name + "_bounds", dimension_name, com, rank);
  dimension.set_string("bounds", dimension_name + "_bounds");

  short_name = name;
  use_bounds = true;
}


//! Read timeseries data from a NetCDF file \c filename.
PetscErrorCode Timeseries::read(const PIO &nc, bool use_reference_date) {
  PetscErrorCode ierr;

  bool exists, found_by_standard_name;
  vector<string> dims;
  string time_name, standard_name = var.get_string("standard_name"),
    name_found;

  ierr = nc.inq_var(short_name, standard_name,
                    exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (!exists) {
    ierr = PetscPrintf(com,
		      "PISM ERROR: Can't find '%s' ('%s') in '%s'.\n",
		       short_name.c_str(), standard_name.c_str(),
                       nc.inq_filename().c_str());
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = nc.inq_vardims(name_found, dims); CHKERRQ(ierr);

  if (dims.size() != 1) {
    ierr = PetscPrintf(com,
		       "PISM ERROR: Variable '%s' in '%s' depends on %d dimensions,\n"
		       "            but a time-series variable can only depend on 1 dimension.\n",
		       short_name.c_str(),
                       nc.inq_filename().c_str(),
                       dims.size()); CHKERRQ(ierr);
    PISMEnd();
  }

  time_name = dims[0];


  dimension.init(time_name, time_name, com, rank);

  ierr = dimension.read(nc, use_reference_date, time); CHKERRQ(ierr);
  bool is_increasing = true;
  for (unsigned int j = 1; j < time.size(); ++j) {
    if (time[j] - time[j-1] < 1e-16) {
      is_increasing = false;
      break;
    }
  }
  if (!is_increasing) {
    ierr = PetscPrintf(com, "PISM ERROR: dimension '%s' has to be strictly increasing (read from '%s').\n",
		       dimension.short_name.c_str(), nc.inq_filename().c_str());
    PISMEnd();
  }

  string time_bounds_name;
  ierr = dimension.get_bounds_name(nc, time_bounds_name); CHKERRQ(ierr);

  if (!time_bounds_name.empty()) {
    use_bounds = true;

    bounds.init(time_bounds_name, time_name, com, rank);
    ierr = bounds.set_units(dimension.get_string("units")); CHKERRQ(ierr);

    ierr = bounds.read(nc, use_reference_date, time_bounds); CHKERRQ(ierr);
  } else {
    use_bounds = false;
  }

  // Do not use the reference date. This may be a problem if someone needs to
  // read a time-series with the meaning of "time depending on time", but this is
  // not likely.
  ierr = var.read(nc, false, values); CHKERRQ(ierr);

  if (time.size() != values.size()) {
    ierr = PetscPrintf(com, "PISM ERROR: variables %s and %s in %s have different numbers of values.\n",
		       dimension.short_name.c_str(),
		       var.short_name.c_str(),
		       nc.inq_filename().c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = var.report_range(values); CHKERRQ(ierr);

  return 0;
}

//! Write timeseries data to a NetCDF file \c filename.
PetscErrorCode Timeseries::write(const PIO &nc) {
  PetscErrorCode ierr;

  // write the dimensional variable; this call should go first
  ierr = dimension.write(nc, 0, time); CHKERRQ(ierr);
  ierr = var.write(nc, 0, values); CHKERRQ(ierr);
  
  if (use_bounds) {
    ierr = bounds.write(nc, 0, time_bounds); CHKERRQ(ierr);
  }

  return 0;
}

//! Get a value of timeseries at time \c t.
/*! Returns the first value or the last value if t is out of range on the left
  and right, respectively.

  Uses time bounds if present (interpreting data as piecewise-constant) and
  uses linear interpolation otherwise.
 */
double Timeseries::operator()(double t) {

  // piecewise-constant case:
  if (use_bounds) {
    vector<double>::iterator j;

    j = lower_bound(time_bounds.begin(), time_bounds.end(), t); // binary search

    if (j == time_bounds.end())
      return values.back(); // out of range (on the right)

    int i = (int)(j - time_bounds.begin());

    if (i == 0)
      return values[0];         // out of range (on the left)

    if (i % 2 == 0) {
      PetscPrintf(com,
                  "PISM ERROR: time bounds array in %s does not represent continguous time intervals.\n"
                  "            (PISM was trying to compute %s at time %3.3f years.)\n",
                  bounds.short_name.c_str(), short_name.c_str(), t);
      PISMEnd();
    }

    return values[(i-1)/2];
  }

  // piecewise-linear case:
  vector<double>::iterator end = time.end(), j;
  
  j = lower_bound(time.begin(), end, t); // binary search

  if (j == end)
    return values.back(); // out of range (on the right)

  int i = (int)(j - time.begin());

  if (i == 0) {
    return values[0];	// out of range (on the left)
  }

  double dt = time[i] - time[i - 1];
  double dv = values[i] - values[i - 1];
  
  return values[i - 1] + (t - time[i - 1]) / dt * dv;
}

//! Get a value of timeseries by index.
/*!
  Stops if the index is out of range.
 */
double Timeseries::operator[](unsigned int j) const {

#if (PISM_DEBUG==1)
  if (j >= values.size()) {
    PetscPrintf(com, "ERROR: Timeseries %s: operator[]: invalid argument: size=%d, index=%d\n",
		var.short_name.c_str(), values.size(), j);
    PISMEnd();
  }
#endif

  return values[j];
}

//! \brief Compute an average of a time-series over interval (t,t+dt) using
//! trapezoidal rule with N sub-intervals.
double Timeseries::average(double t, double dt, unsigned int N) {
  vector<double> V(N+1);
 
  for (unsigned int i = 0; i < N+1; ++i) {
    double t_i = t + (dt / N) * i;
    V[i] = (*this)(t_i);
  }

  double sum = 0;
  for (unsigned int i = 0; i < N; ++i) {
    sum += V[i] + V[i+1];
  }

  return sum / (2*N);
}

//! Append a pair (t,v) to the timeseries.
PetscErrorCode Timeseries::append(double v, double a, double b) {
  time.push_back(b);
  values.push_back(v);
  time_bounds.push_back(a);
  time_bounds.push_back(b);
  return 0;
}


//! Set the internal units for the values of a time-series.
PetscErrorCode Timeseries::set_units(string units, string glaciological_units) {
  if (!units.empty())
    var.set_units(units);
  if (!glaciological_units.empty())
    var.set_glaciological_units(glaciological_units);
  return 0;
}

//! Set the internal units for the dimension variable of a time-series.
PetscErrorCode Timeseries::set_dimension_units(string units, string glaciological_units) {
  if (!units.empty())
    dimension.set_units(units);
  if (!glaciological_units.empty())
    dimension.set_glaciological_units(glaciological_units);
  return 0;
}

//! Set a string attribute.
PetscErrorCode Timeseries::set_attr(string name, string value) {
  var.set_string(name, value);
  return 0;
}

//! Get a string attribute.
string Timeseries::get_string(string name) {
  return var.get_string(name);
}

//! Set a single-valued scalar attribute.
PetscErrorCode Timeseries::set_attr(string name, double value) {
  var.set(name, value);
  return 0;
}

//! Returns the length of the time-series stored.
/*!
  This length is changed by read() and append().
 */
int Timeseries::length() {
  return (int)values.size();
}


//----- DiagnosticTimeseries

DiagnosticTimeseries::DiagnosticTimeseries(IceGrid *g, string name, string dimension_name)
  : Timeseries(g, name, dimension_name) {

  buffer_size = (size_t)g->config.get("timeseries_buffer_size");
  start = 0;
  rate_of_change = false;
  dimension.set_string("calendar", g->config.get_string("calendar"));
  dimension.set_string("long_name", "time");
  dimension.set_string("axis", "T");
}

//! Destructor; makes sure that everything is written to a file.
DiagnosticTimeseries::~DiagnosticTimeseries() {
  flush();
}

//! Adds the (t,v) pair to the interpolation buffer.
/*! The interpolation buffer holds 2 values only (for linear interpolation).
 *
 * If this DiagnosticTimeseries object is reporting a "rate of change",
 * append() has to be called with the "cumulative" quantity as the V argument.
 */
PetscErrorCode DiagnosticTimeseries::append(double V, double /*a*/, double b) {

  if (rate_of_change && v.empty()) {
    v_previous = V;
  }

  // append to the interpolation buffer
  t.push_back(b);
  v.push_back(V);

  if (t.size() == 3) {
    t.pop_front();
    v.pop_front();
  }

  return 0;
}

//! \brief Use linear interpolation to find the value of a scalar diagnostic
//! quantity at time \c T and store the obtained pair (T, value).
PetscErrorCode DiagnosticTimeseries::interp(double a, double b) {
  PetscErrorCode ierr;

  if (t.empty()) {
    SETERRQ(com, 1, "DiagnosticTimeseries::interp(...): interpolation buffer is empty");
  }

  if (t.size() == 1) {
    time.push_back(b);
    values.push_back(GSL_NAN);
    time_bounds.push_back(a);
    time_bounds.push_back(b);
    return 0;
  }

  if ((b < t[0]) || (b > t[1])) {
    SETERRQ1(com, 1, "DiagnosticTimeseries::interp(...): requested time %f is not within the last time-step!",
             b);
  }

  double 
    // compute the "cumulative" quantity using linear interpolation
    v_current = v[0] + (b - t[0]) / (t[1] - t[0]) * (v[1] - v[0]),
    // the value to report
    value = v_current;

  if (rate_of_change) {
    // use backward-in-time finite difference to compute the rate of change:
    value = (v_current - v_previous) / (b - a);

    // remember the value of the "cumulative" quantity for differencing during
    // the next call:
    v_previous = v_current;
  }

  // use the right endpoint as the 'time' record (the midpoint is also an option)
  time.push_back(b);
  values.push_back(value);

  // save the time bounds
  time_bounds.push_back(a);
  time_bounds.push_back(b);

  if (time.size() == buffer_size) {
    ierr = flush(); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode DiagnosticTimeseries::init(string filename) {
  PetscErrorCode ierr;
  PIO nc(com, rank, "netcdf3"); // OK to use netcdf3
  unsigned int len = 0;

  // Get the number of records in the file (for appending):
  int file_exists = 0;
  if (rank == 0) {
    // Check if the file exists:
    if (FILE *f = fopen(filename.c_str(), "r")) {
      file_exists = 1;
      fclose(f);
    } else {
      file_exists = 0;
    }
  }
  ierr = MPI_Bcast(&file_exists, 1, MPI_INT, 0, com); CHKERRQ(ierr);

  if (file_exists == 1) {
    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_dimlen(dimension.short_name, len); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
  }

  output_filename = filename;
  start = len;

  return 0;
}


  //! Writes data to a file.
PetscErrorCode DiagnosticTimeseries::flush() {
  PetscErrorCode ierr;
  PIO nc(com, rank, "netcdf3"); // OK to use netcdf3
  unsigned int len = 0;

  // return cleanly if this DiagnosticTimeseries object was created but never
  // used:
  if (output_filename.empty())
    return 0;

  if (time.empty())
    return 0;

  ierr = nc.open(output_filename, PISM_WRITE, true); CHKERRQ(ierr);
  ierr = nc.inq_dimlen(dimension.short_name, len); CHKERRQ(ierr);

  if (len > 0) {
    double last_time;
    ierr = nc.inq_dim_limits(dimension.dimension_name, NULL, &last_time); CHKERRQ(ierr);
    if (last_time < time.front()) {
      start = len;
    }
  }

  if (len == (unsigned int)start) {
    ierr = dimension.write(nc, start, time);   CHKERRQ(ierr);
    ierr = bounds.write(nc, start, time_bounds);   CHKERRQ(ierr);
  }
  ierr = var.write(nc, start, values); CHKERRQ(ierr);

  start += time.size();

  time.clear();
  values.clear();
  time_bounds.clear();

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

void DiagnosticTimeseries::reset() {
  time.clear();
  values.clear();
  time_bounds.clear();
  start = 0;
  t.clear();
  v.clear();
}

