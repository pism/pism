// Copyright (C) 2011 Constantine Khroulev
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

  //! \brief Intialize using command-line options.
  virtual PetscErrorCode init();

  //! \brief Sets the current time (in seconds since the reference time).
  virtual void set(PetscReal new_time)
  { time_in_seconds = new_time; }

  virtual void set_start(PetscReal new_start)
  { run_start = new_start; }

  virtual void set_end(PetscReal new_end)
  { run_end = new_end; }

  //! \brief Sets the reference date string.
  virtual void set_reference_date(string str)
  { reference_date = str; }

  //! \brief Advance by delta_t seconds.
  virtual void step(PetscReal delta_t) 
  { time_in_seconds += delta_t; }
  
  //! \brief Current time, in seconds.
  virtual PetscReal current()
  { return time_in_seconds; }

  virtual PetscReal start()
  { return run_start; }

  virtual PetscReal end()
  { return run_end; }

  //! \brief Returns time since the origin modulo period.
  virtual PetscReal mod(PetscReal time, PetscReal period)
  {
    if (period <= 0)
      return time;

    return time - floor(time / period) * period;
  }

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year.
  virtual PetscReal year_fraction(PetscReal T)
  { return seconds_to_years(T) - floor(seconds_to_years(T)); }

  //! \brief Returns the year corresponding to time T.
  virtual PetscReal year(PetscReal T)
  { return seconds_to_years(T); }

  //! \brief Returns current time, in years.
  virtual PetscReal year()
  { return time_in_seconds / secpera; }

  //! \brief All times are interpreted as "time since the reference time", even
  //! if this implementation does not take advantage of that.
  virtual PetscReal start_year()
  { return run_start / secpera; }

  //! \brief All times are interpreted as "time since the reference time", even
  //! if this implementation does not take advantage of that.
  virtual PetscReal end_year()
  { return run_end / secpera; }

  //! \brief Returns the length of the current run, in years.
  virtual PetscReal run_length_years()
  { return (run_end - run_start) / secpera; }

  //! \brief Returns the CF- (and UDUNITS) compliant units string.
  virtual string units()
  { return string("seconds since ") + reference_date; }

  //! \brief Returns the calendar string.
  virtual string calendar()
  { return calendar_string; }

  virtual PetscReal seconds_to_years(PetscReal T)
  { return T / secpera; }

  virtual PetscReal years_to_seconds(PetscReal T)
  { return T * secpera; }

protected:
  MPI_Comm com;
  const NCConfigVariable &config;
  PetscReal secpera;      //!< number of seconds in a year, for unit conversion
  PetscReal time_in_seconds, //!< current time, in seconds since the reference time
    run_start,                  //!< run start time, in seconds since the reference time
    run_end;                    //!< run end tim, in seconds since the reference time
  string reference_date,     //!< CF reference date; used in the units string
    calendar_string;         //!< CF calendard string
};

#endif /* _PISMTIME_H_ */
