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

#endif // __Timeseries_hh
