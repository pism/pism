// Copyright (C) 2009--2011 Constantine Khroulev
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

Timeseries::Timeseries(IceGrid *g, string name, string dimension_name) {
  com = g->com;
  rank = g->rank;

  dimension.init(dimension_name, dimension_name, com, rank);
  var.init(name, dimension_name, com, rank);

  short_name = name;
}

Timeseries::Timeseries(MPI_Comm c, PetscMPIInt r, string name, string dimension_name) {
  com = c;
  rank = r;
  dimension.init(dimension_name, dimension_name, com, rank);
  var.init(name, dimension_name, com, rank);

  short_name = name;
}

//! Read timeseries data from a NetCDF file \c filename.
PetscErrorCode Timeseries::read(const char filename[]) {
  PetscErrorCode ierr;

  NCTool nc(com, rank);
  bool exists;
  vector<int> dimids;
  int varid;
  string time_name, standard_name = var.get_string("standard_name");
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  ierr = nc.find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);

  if (!exists) {
    ierr = PetscPrintf(com,
		      "PISM ERROR: Can't find '%s' ('%s') in '%s'.\n",
		       short_name.c_str(), standard_name.c_str(), filename);
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = nc.inq_dimids(varid, dimids); CHKERRQ(ierr);

  if (dimids.size() != 1) {
    ierr = PetscPrintf(com,
		       "PISM ERROR: Variable '%s' in '%s' depends on %d dimensions,\n"
		       "            but a time-series variable can only depend on 1 dimension.\n",
		       short_name.c_str(), filename, dimids.size()); CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = nc.inq_dimname(dimids[0], time_name); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  dimension.init(time_name, time_name, com, rank);

  ierr = dimension.read(filename, time); CHKERRQ(ierr);
  bool is_increasing = true;
  for (unsigned int j = 1; j < time.size(); ++j) {
    if (time[j] - time[j-1] < 1e-16) {
      is_increasing = false;
      break;
    }
  }
  if (!is_increasing) {
    ierr = PetscPrintf(com, "PISM ERROR: dimension '%s' has to be strictly increasing (read from '%s').\n",
		       dimension.short_name.c_str(), filename);
    PISMEnd();
  }

  ierr = var.read(filename, values); CHKERRQ(ierr);

  if (time.size() != values.size()) {
    ierr = PetscPrintf(com, "PISM ERROR: variables %s and %s in %s have different numbers of values.\n",
		       dimension.short_name.c_str(),
		       var.short_name.c_str(),
		       filename); CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = var.report_range(values); CHKERRQ(ierr);

  return 0;
}

//! Write timeseries data to a NetCDF file \c filename.
PetscErrorCode Timeseries::write(const char filename[]) {
  PetscErrorCode ierr;

  // write the dimensional variable; this call should go first
  ierr = dimension.write(filename, 0, time); CHKERRQ(ierr);
  ierr = var.write(filename, 0, values); CHKERRQ(ierr);

  return 0;
}

//! Get an (linearly interpolated) value of timeseries at time \c t.
/*! Returns the first value or the last value if t is out of range on the left
  and right, respectively.
 */
double Timeseries::operator()(double t) {

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

  if (j >= values.size()) {
    PetscPrintf(com, "ERROR: Timeseries %s: operator[]: invalid argument: size=%d, index=%d\n",
		var.short_name.c_str(), values.size(), j);
    PISMEnd();
  }

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
PetscErrorCode Timeseries::append(double t, double v) {
  time.push_back(t);
  values.push_back(v);
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

  buffer_size = 10000;		// just a default
  start = 0;
}

DiagnosticTimeseries::DiagnosticTimeseries(MPI_Comm c, PetscMPIInt r, string name, string dimension_name)
  : Timeseries(c, r, name, dimension_name) {

  buffer_size = 10000;		// just a default
  start = 0;
}

//! Destructor; makes sure that everything is written to a file.
DiagnosticTimeseries::~DiagnosticTimeseries() {
  flush();
}

//! Adds the (t,v) pair to the interpolation buffer.
/*! The interpolation buffer holds 2 values only (for linear interpolation).
 */
PetscErrorCode DiagnosticTimeseries::append(double T, double V) {
  // append to the interpolation buffer:
  t.push_back(T);
  v.push_back(V);

  if (t.size() == 3) {
    t.pop_front();
    v.pop_front();
  }

  return 0;
}

//! \brief Use linear interpolation to find the value of a scalar diagnostic
//! quantity at time \c T and store the obtained pair (T, value).
PetscErrorCode DiagnosticTimeseries::interp(double T) {
  PetscErrorCode ierr;

  if (t.empty()) {
    SETERRQ(1, "DiagnosticTimeseries::interp(...): interpolation buffer is empty");
  }

  // Special case: allow requests that match the (only) available time.
  if (t.size() == 1) {
    if (PetscAbs(T - t[0]) < 1e-8) {
      time.push_back(T);
      values.push_back(v[0]);
    }
    return 0;
  }

  if ((T < t[0]) || (T > t[1])) {
    SETERRQ1(1, "DiagnosticTimeseries::interp(...): requested time %f is not within the last time-step!",
	     T);
  }

  // linear interpolation:
  double tmp = v[0] + (T - t[0]) / (t[1] - t[0]) * (v[1] - v[0]);

  time.push_back(T);
  values.push_back(tmp);

  if (time.size() == buffer_size) {
    ierr = flush(); CHKERRQ(ierr);
  }

  return 0;
}

//! Writes data to a file.
PetscErrorCode DiagnosticTimeseries::flush() {
  PetscErrorCode ierr;
  NCTool nc(com, rank);
  unsigned int len = 0;

  // return cleanly if this DiagnosticTimeseries object was created but never
  // used:
  if (output_filename.empty())
    return 0;

  // Get the length of the time dimension; if other time-series object was
  // written at this time-step already, len will be greater than start, so we
  // don't need to write the dimensional variable.
  ierr = nc.open_for_reading(output_filename.c_str()); CHKERRQ(ierr);
  ierr = nc.get_dim_length(dimension.short_name.c_str(), &len); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  if (len > 0) {
    double last_time;
    ierr = nc.get_dim_limits(dimension.dimension_name, NULL, &last_time); CHKERRQ(ierr);
    if (last_time/secpera < time.front()) {
      start = len;
    }
  }

  if (len == (unsigned int)start) {
    ierr = dimension.write(output_filename.c_str(), start, time);   CHKERRQ(ierr);
  }
  ierr = var.write(output_filename.c_str(), start, values); CHKERRQ(ierr);

  start += time.size();
  time.clear();
  values.clear();

  return 0;
}
