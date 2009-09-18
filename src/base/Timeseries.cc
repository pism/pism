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

#include "Timeseries.hh"
#include <algorithm>

Timeseries::Timeseries(IceGrid *g, string name, string dimension_name) {
  com = g->com;
  rank = g->rank;

  dimension.init(dimension_name, com, rank);
  dimension.dimension_name = dimension_name;

  var.init(name, com, rank);
  var.dimension_name = dimension_name;
}

Timeseries::Timeseries(MPI_Comm c, PetscMPIInt r, string name, string dimension_name) {
  com = c;
  rank = r;
  dimension.init(dimension_name, com, rank);
  var.init(name, com, rank);
}

//! Read timeseries data from a NetCDF file \c filename.
PetscErrorCode Timeseries::read(const char filename[]) {
  PetscErrorCode ierr;

  ierr = dimension.read(filename, time); CHKERRQ(ierr);
  ierr = var.read(filename, values); CHKERRQ(ierr);

  if (time.size() != values.size()) {
    ierr = PetscPrintf(com, "PISM ERROR: variables %s and %s in %s have different numbers of values.\n",
		       dimension.short_name.c_str(),
		       var.short_name.c_str(),
		       filename); CHKERRQ(ierr);
    PetscEnd();
  }

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

  int i = j - time.begin();

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
    PetscEnd();
  }

  return values[j];
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
  return values.size();
}


//----- DiagnosticTimeseries

DiagnosticTimeseries::DiagnosticTimeseries(IceGrid *g, string name, string dimension_name)
  : Timeseries(g, name, dimension_name) {

  buffer_size = 10000;		// just a default
  start = 0;
  output_filename = name + ".nc";
}

DiagnosticTimeseries::DiagnosticTimeseries(MPI_Comm c, PetscMPIInt r, string name, string dimension_name)
  : Timeseries(c, r, name, dimension_name) {

  buffer_size = 10000;		// just a default
  start = 0;
  output_filename = name + ".nc";
}

DiagnosticTimeseries::~DiagnosticTimeseries() {
  flush();
}

PetscErrorCode DiagnosticTimeseries::append(double t, double v) {
  PetscErrorCode ierr;

  time.push_back(t);
  values.push_back(v);

  if (time.size() == buffer_size) {
    ierr = flush(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode DiagnosticTimeseries::set_output_prefix(string prefix) {
  PetscErrorCode ierr;
  NCTool nc(com, rank);

  output_filename = prefix + var.short_name + ".nc";

  // This will move the file aside if it exists already.
  ierr = nc.open_for_writing(output_filename.c_str(), false, false); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode DiagnosticTimeseries::flush() {
  PetscErrorCode ierr;

  // write the dimensional variable; this call should go first
  ierr = dimension.write(output_filename.c_str(), start, time);   CHKERRQ(ierr);
  ierr =       var.write(output_filename.c_str(), start, values); CHKERRQ(ierr);

  start += buffer_size;
  time.clear();
  values.clear();

  return 0;
}
