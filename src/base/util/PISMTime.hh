// Copyright (C) 2011, 2012 Constantine Khroulev
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
  PISMTime(MPI_Comm c, const NCConfigVariable &conf);
  virtual ~PISMTime() {}

  //! \brief Sets the current time (in seconds since the reference time).
  void set(double new_time)
  { time_in_seconds = new_time; }

  void set_start(double new_start)
  { run_start = new_start; }

  void set_end(double new_end)
  { run_end = new_end; }

  //! \brief Sets the reference date string.
  void set_reference_date(string str)
  { reference_date = str; }

  //! \brief Advance by delta_t seconds.
  void step(double delta_t)
  {
    time_in_seconds += delta_t;

    // If we are less than 0.001 second from the end of the run, reset
    // time_in_seconds to avoid taking a very small (and useless) time step.
    if (run_end > time_in_seconds &&
        run_end - time_in_seconds < 1e-3)
      time_in_seconds = run_end;
  }

  //! \brief Current time, in seconds.
  double current()
  { return time_in_seconds; }

  double start()
  { return run_start; }

  double end()
  { return run_end; }

  //! \brief Returns the CF- (and UDUNITS) compliant units string.
  /*!
   * This units string is saved in the output file. Always contains a reference
   * date, even if it is not used by PISM.
   */
  string CF_units()
  { return string("seconds since ") + reference_date; }

  //! \brief Returns the calendar string.
  string calendar()
  { return calendar_string; }

  //! \brief Returns the length of the current run, in years.
  string run_length();

  double seconds_to_years(double T)
  { return T / secpera; }

  double years_to_seconds(double T)
  { return T * secpera; }

  // Virtual methods:

  //! \brief Intialize using command-line options.
  virtual PetscErrorCode init();

  //! \brief Internal time units.
  /*!
   * May or may not contain a reference date. (The base class PISMTime does not
   * use the reference date, while PISMGregorianTime does.)
   */
  virtual string units()
  { return "seconds"; }

  virtual bool use_reference_date()
  { return false; }

  //! \brief Returns time since the origin modulo period.
  virtual double mod(double time, double period);

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year. Only useful in codes with a "yearly cycle" (such as the PDD model).
  virtual double year_fraction(double T);

  //! \brief Returns the date corresponding to time T.
  virtual string date(double T);

  //! \brief Returns current time, in years. Only for reporting.
  virtual string date();

  //! \brief All times are interpreted as "time since the reference time", even
  //! if this implementation does not take advantage of that.
  virtual string start_date();

  //! \brief All times are interpreted as "time since the reference time", even
  //! if this implementation does not take advantage of that.
  virtual string end_date();

protected:
  MPI_Comm com;
  const NCConfigVariable &config;
  double secpera;      //!< number of seconds in a year, for unit conversion
  double time_in_seconds, //!< current time, in seconds since the reference time
    run_start,                  //!< run start time, in seconds since the reference time
    run_end;                    //!< run end tim, in seconds since the reference time
  string reference_date,     //!< CF reference date; used in the units string
    calendar_string;         //!< CF calendard string
};

#endif /* _PISMTIME_H_ */
