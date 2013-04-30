// Copyright (C) 2011, 2012, 2013 Constantine Khroulev
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

#ifndef _PISMTIME_H_
#define _PISMTIME_H_

#include "pism_const.hh"
#include "NCVariable.hh"

//! \brief Time management class.
/*!
 * This is to make it possible to switch between different implementations.
 *
 * For example: 365-day no-leap calendar for spinups
 * Gregorian calendar for XX-century forcing runs
 *
 * This base class implements the 365.24-day no-leap version.
 *
 * We want to be able to count time since a particular date, so it is helpful
 * to keep in mind that the year "1986" in this context is not the year of the
 * Chernobyl disaster but a year 1986 years since some date.
 */
class PISMTime
{
public:
  PISMTime(MPI_Comm c, const NCConfigVariable &conf, PISMUnitSystem units_system);
  virtual ~PISMTime();

  //! \brief Sets the current time (in seconds since the reference time).
  void set(double new_time);

  void set_start(double new_start);

  void set_end(double new_end);

  //! \brief Advance by delta_t seconds.
  void step(double delta_t);

  //! \brief Current time, in seconds.
  double current();

  double start();

  double end();

  //! \brief Returns the calendar string.
  string calendar();

  //! \brief Returns the length of the current run, in years.
  string run_length();

  double seconds_to_years(double T);
  double years_to_seconds(double T);

  // Virtual methods:

  //! \brief Intialize using command-line options.
  virtual PetscErrorCode init();

  PetscErrorCode parse_times(string spec, vector<double> &result);


  //! \brief Returns the CF- (and UDUNITS) compliant units string.
  /*!
   * This units string is saved in the output file. Always contains a reference
   * date, even if it is not used by PISM.
   */
  virtual string CF_units_string();

  //! \brief Internal time units.
  /*!
   * May or may not contain a reference date. (The base class PISMTime does not
   * use the reference date, while PISMTime_Calendar does.)
   */
  virtual string units_string();

  virtual bool use_reference_date();

  //! \brief Returns time since the origin modulo period.
  virtual double mod(double time, double period);

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year. Only useful in codes with a "yearly cycle" (such as the PDD model).
  virtual double year_fraction(double T);

  //! \brief Returns the model time in seconds corresponding to the
  //! beginning of the year `T` falls into.
  virtual double calendar_year_start(double T);

  //! Increment time `T` by a given amount and return resulting model
  //! time in seconds.
  virtual double increment_date(double T, int years);

  //! \brief Returns the date corresponding to time T.
  virtual string date(double T);

  //! \brief Returns current time, in years. Only for reporting.
  virtual string date();

  //! Date corresponding to the beginning of the run.
  virtual string start_date();

  //! Date corresponding to the end of the run.
  virtual string end_date();

protected:
  PetscErrorCode parse_list(string spec, vector<double> &result);

  virtual PetscErrorCode process_ys(double &result, bool &flag);
  virtual PetscErrorCode process_y(double &result, bool &flag);
  virtual PetscErrorCode process_ye(double &result, bool &flag);

  virtual PetscErrorCode compute_times(double time_start, double delta, double time_end,
                                       string keyword,
                                       vector<double> &result);

  PetscErrorCode compute_times_simple(double time_start, double delta, double time_end,
                                      vector<double> &result);

  virtual PetscErrorCode parse_range(string spec, vector<double> &result);

  virtual PetscErrorCode parse_date(string spec, double *result);

  virtual PetscErrorCode parse_interval_length(string spec, string &keyword, double *result);

  MPI_Comm m_com;
  const NCConfigVariable &m_config;
  PISMUnitSystem m_unit_system;
  PISMUnit m_time_units;
  double m_secpera;      //!< number of seconds in a year, for unit conversion
  double m_time_in_seconds, //!< current time, in seconds since the reference time
    m_run_start,                  //!< run start time, in seconds since the reference time
    m_run_end;                    //!< run end tim, in seconds since the reference time
  string m_calendar_string;       //!< CF calendar string
};

#endif /* _PISMTIME_H_ */
