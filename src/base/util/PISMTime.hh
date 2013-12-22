// Copyright (C) 2011, 2012, 2013, 2014 Constantine Khroulev
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

/**
 * Returns 0 if `name` is a name of a supported calendar, 1 otherwise.
 */
inline bool pism_is_valid_calendar_name(std::string name) {
  // Calendar names from the CF Conventions document (except the
  // 366_day (all_leap)):
  if (name == "standard"            ||
      name == "gregorian"           ||
      name == "proleptic_gregorian" ||
      name == "noleap"              ||
      name == "365_day"             ||
      name == "julian"              ||
      name == "360_day"             ||
      name == "none") {
    return true;
  }

  return false;
}

//! \brief Time management class.
/*!
 * This is to make it possible to switch between different implementations.
 *
 * For example: 365-day no-leap calendar for spinups
 * Gregorian calendar for XX-century forcing runs
 *
 * This base class implements the 365-day no-leap version.
 *
 * We want to be able to count time since a particular date, so it is helpful
 * to keep in mind that the year "1986" in this context is not the year of the
 * Chernobyl disaster but a year 1986 years since some date.
 */
class PISMTime
{
public:
  PISMTime(MPI_Comm c, const PISMConfig &conf, std::string calendar,
           PISMUnitSystem units_system);
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
  std::string calendar();

  //! \brief Returns the length of the current run, in years.
  std::string run_length();

  // Virtual methods:

  //! \brief Intialize using command-line options.
  virtual PetscErrorCode init();

  void init_calendar(std::string calendar);

  PetscErrorCode parse_times(std::string spec, std::vector<double> &result);

  //! \brief Returns the CF- (and UDUNITS) compliant units string.
  /*!
   * This units string is saved in the output file. Always contains a reference
   * date, even if it is not used by PISM.
   */
  virtual std::string CF_units_string();

  //! \brief Internal time units.
  /*!
   * May or may not contain a reference date. (The base class PISMTime does not
   * use the reference date, while PISMTime_Calendar does.)
   */
  virtual std::string units_string();

  virtual std::string CF_units_to_PISM_units(std::string input);

  //! \brief Returns time since the origin modulo period.
  virtual double mod(double time, unsigned int period_years);

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year. Only useful in codes with a "yearly cycle" (such as the PDD model).
  virtual double year_fraction(double T);

  //! \brief Convert the day number to the year fraction.
  virtual double day_of_the_year_to_day_fraction(unsigned int day);

  //! \brief Returns the model time in seconds corresponding to the
  //! beginning of the year `T` falls into.
  virtual double calendar_year_start(double T);

  //! Increment time `T` by a given amount and return resulting model
  //! time in seconds.
  virtual double increment_date(double T, int years);

  //! \brief Returns the date corresponding to time T.
  virtual std::string date(double T);

  //! \brief Returns current time, in years. Only for reporting.
  virtual std::string date();

#if (PISM_DEBUG==1)
  //! \brief Returns current time, in years. Only for debugging.
  virtual double current_years() {
    return seconds_to_years(current());
  }
#endif

  //! Date corresponding to the beginning of the run.
  virtual std::string start_date();

  //! Date corresponding to the end of the run.
  virtual std::string end_date();

  //! @brief Convert time interval from seconds to given units. Handle
  //! 'years' using the year length corresponding to the calendar.
  virtual double convert_time_interval(double T, std::string units);

protected:
  PetscErrorCode parse_list(std::string spec, std::vector<double> &result);

  virtual PetscErrorCode process_ys(double &result, bool &flag);
  virtual PetscErrorCode process_y(double &result, bool &flag);
  virtual PetscErrorCode process_ye(double &result, bool &flag);

  virtual PetscErrorCode compute_times(double time_start, double delta, double time_end,
                                       std::string keyword,
                                       std::vector<double> &result);

  PetscErrorCode compute_times_simple(double time_start, double delta, double time_end,
                                      std::vector<double> &result);

  virtual PetscErrorCode parse_range(std::string spec, std::vector<double> &result);

  virtual PetscErrorCode parse_date(std::string spec, double *result);

  virtual PetscErrorCode parse_interval_length(std::string spec, std::string &keyword,
                                               double *result);

  double years_to_seconds(double input);
  double seconds_to_years(double input);

protected:
  MPI_Comm m_com;
  const PISMConfig &m_config;
  PISMUnitSystem m_unit_system;
  PISMUnit m_time_units;
  double m_year_length;      //!< number of seconds in a year, for "mod" and "year fraction"
  double m_time_in_seconds, //!< current time, in seconds since the reference time
    m_run_start,                  //!< run start time, in seconds since the reference time
    m_run_end;                    //!< run end tim, in seconds since the reference time
  std::string m_calendar_string;       //!< CF calendar string
};

#endif /* _PISMTIME_H_ */
