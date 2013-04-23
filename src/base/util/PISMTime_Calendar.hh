// Copyright (C) 2012, 2013 PISM Authors
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

#ifndef _PISMGREGORIANTIME_H_
#define _PISMGREGORIANTIME_H_

#include "PISMTime.hh"
#include "PISMUnits.hh"

class PISMTime_Calendar : public PISMTime
{
public:
  PISMTime_Calendar(MPI_Comm c, const NCConfigVariable &conf,
                    PISMUnitSystem units_system);
  virtual ~PISMTime_Calendar();

  virtual PetscErrorCode init();

  virtual PetscErrorCode init_from_file(string filename);

  virtual double mod(double time, double period);

  virtual double year_fraction(double T);

  virtual string date(double T);

  virtual string date();

  virtual string start_date();

  virtual string end_date();

  virtual string units_string()
  { return CF_units_string(); }

  virtual bool use_reference_date()
  { return true; }

  virtual double calendar_year_start(double T);

  virtual double increment_date(double T, int years);

protected:
  virtual PetscErrorCode compute_times(double time_start, double delta, double time_end,
                                       string keyword,
                                       vector<double> &result);

  virtual PetscErrorCode process_ys(double &result, bool &flag);
  virtual PetscErrorCode process_y(double &result, bool &flag);
  virtual PetscErrorCode process_ye(double &result, bool &flag);

  virtual PetscErrorCode parse_date(string spec, double *result);

  virtual PetscErrorCode parse_interval_length(string spec, string &keyword, double *result);

  PetscErrorCode compute_times_monthly(vector<double> &result);

  PetscErrorCode compute_times_yearly(vector<double> &result);
private:
  // Hide copy constructor / assignment operator.
  PISMTime_Calendar(PISMTime_Calendar const &);
  PISMTime_Calendar & operator=(PISMTime_Calendar const &);
};


#endif /* _PISMGREGORIANTIME_H_ */
