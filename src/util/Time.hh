// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 Constantine Khroulev
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

#include <vector>
#include <memory>

#include "pism/util/pism_utilities.hh"
#include "pism/util/Units.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {

/**
 * Returns 0 if `name` is a name of a supported calendar, 1 otherwise.
 */
inline bool pism_is_valid_calendar_name(const std::string &name) {
  // Calendar names from the CF Conventions document (except the
  // 366_day (all_leap)):
  return member(name, {"standard", "gregorian", "proleptic_gregorian",
                       "noleap", "365_day", "julian", "360_day"});
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
class Time
{
public:
  Time(MPI_Comm com, Config::ConstPtr config,
       const Logger &log,
       units::System::Ptr unit_system);
  virtual ~Time() = default;

  typedef std::shared_ptr<Time> Ptr;
  typedef std::shared_ptr<const Time> ConstPtr;

  //! \brief Sets the current time (in seconds since the reference time).
  void set(double new_time);

  void set_start(double new_start);

  void set_end(double new_end);

  //! \brief Advance by delta_t seconds.
  void step(double delta_t);

  //! \brief Current time, in seconds.
  double current() const;

  double start() const;

  double end() const;

  //! \brief Returns the calendar string.
  std::string calendar() const;

  //! \brief Returns the length of the current run, in years.
  std::string run_length() const;

  void init_calendar(const std::string &calendar);

  std::vector<double> parse_times(const std::string &spec) const;

  //! \brief Internal time units.
  /*!
   * May or may not contain a reference date. (The base class Time does not
   * use the reference date, while Time_Calendar does.)
   */
  std::string units_string() const;

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year. Only useful in codes with a "yearly cycle" (such as the PDD model).
  double year_fraction(double T) const;

  //! \brief Convert the day number to the year fraction.
  double day_of_the_year_to_day_fraction(unsigned int day) const;

  //! \brief Returns the model time in seconds corresponding to the
  //! beginning of the year `T` falls into.
  double calendar_year_start(double T) const;

  //! Increment time `T` by a given amount and return resulting model
  //! time in seconds.
  double increment_date(double T, double years) const;

  //! \brief Returns the date corresponding to time T.
  std::string date(double T) const;

  //! @brief Convert time interval from seconds to given units. Handle
  //! 'years' using the year length corresponding to the calendar.
  double convert_time_interval(double T, const std::string &units) const;

  //! Convert time interval length in years into seconds using the year length
  //! corresponding to the chosen calendar.
  double years_to_seconds(double input) const;

  //! Convert time interval length in seconds into years using the year length
  //! corresponding to the chosen calendar.
  double seconds_to_years(double input) const;
protected:

  std::vector<double> parse_list(const std::string &spec) const;
  std::vector<double> parse_range(const std::string &spec) const;

  void compute_times_simple(double time_start, double delta, double time_end,
                            std::vector<double> &result) const;

  enum IntervalType {YEARLY, MONTHLY, SIMPLE};

  struct Interval {
    double dt;
    IntervalType type;
  };

  void compute_times(double time_start, double time_end,
                     const Interval &interval,
                     std::vector<double> &result) const;

  Interval parse_interval_length(const std::string &spec) const;

protected:
  const Config::ConstPtr m_config;
  const units::System::Ptr m_unit_system;
  units::Unit m_time_units;
  //! number of seconds in a year, for "mod" and "year fraction"
  double m_year_length;

  //! current time, in seconds since the reference time
  double m_time_in_seconds;

  //! run start time, in seconds since the reference time
  double m_run_start;

  //! run end tim, in seconds since the reference time
  double m_run_end;

  //! CF calendar string
  std::string m_calendar_string;
  // True if the calendar has constant year lengths (360_day, 365_day)
  bool m_simple_calendar;

  void init_from_file(MPI_Comm com, const std::string &filename, const Logger &log,
                      bool set_start_time);

  void compute_times_monthly(std::vector<double> &result) const;

  void compute_times_yearly(std::vector<double> &result) const;
};

void check_forcing_duration(const Time &time,
                            double forcing_start,
                            double forcing_end);


} // end of namespace pism

#endif /* _PISMTIME_H_ */
