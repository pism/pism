// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev
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

#include "pism_utilities.hh"
#include "Units.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {

class Config;

std::string calendar_from_options(MPI_Comm com, const Config& config);

/**
 * Returns 0 if `name` is a name of a supported calendar, 1 otherwise.
 */
inline bool pism_is_valid_calendar_name(const std::string &name) {
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
class Time
{
public:
  Time(Config::ConstPtr conf,
       const std::string &calendar,
       units::System::Ptr units_system);
  virtual ~Time();

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

  // Virtual methods:

  //! \brief Intialize using command-line options.
  virtual void init(const Logger &log);

  virtual void init_from_input_file(const PIO &nc,
                                    const std::string &time_name,
                                    const Logger &log);

  void init_calendar(const std::string &calendar);

  std::vector<double> parse_times(const std::string &spec) const;

  //! \brief Returns the CF- (and UDUNITS) compliant units string.
  /*!
   * This units string is saved in the output file. Always contains a reference
   * date, even if it is not used by PISM.
   */
  virtual std::string CF_units_string() const;

  //! \brief Internal time units.
  /*!
   * May or may not contain a reference date. (The base class Time does not
   * use the reference date, while Time_Calendar does.)
   */
  virtual std::string units_string() const;

  virtual std::string CF_units_to_PISM_units(const std::string &input) const;

  //! \brief Returns time since the origin modulo period.
  virtual double mod(double time, unsigned int period_years) const;

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year. Only useful in codes with a "yearly cycle" (such as the PDD model).
  virtual double year_fraction(double T) const;

  //! \brief Convert the day number to the year fraction.
  virtual double day_of_the_year_to_day_fraction(unsigned int day) const;

  //! \brief Returns the model time in seconds corresponding to the
  //! beginning of the year `T` falls into.
  virtual double calendar_year_start(double T) const;

  //! Increment time `T` by a given amount and return resulting model
  //! time in seconds.
  virtual double increment_date(double T, int years) const;

  //! \brief Returns the date corresponding to time T.
  virtual std::string date(double T) const;

  //! \brief Returns current time, in years. Only for reporting.
  virtual std::string date() const;

#if (PISM_DEBUG==1)
  //! \brief Returns current time, in years. Only for debugging.
  virtual double current_years() const {
    return seconds_to_years(current());
  }
#endif

  //! Date corresponding to the beginning of the run.
  virtual std::string start_date() const;

  //! Date corresponding to the end of the run.
  virtual std::string end_date() const;

  //! @brief Convert time interval from seconds to given units. Handle
  //! 'years' using the year length corresponding to the calendar.
  virtual double convert_time_interval(double T, const std::string &units) const;

protected:
  double years_to_seconds(double input) const;
  double seconds_to_years(double input) const;

  std::vector<double> parse_list(const std::string &spec) const;
  std::vector<double> parse_range(const std::string &spec) const;

  void compute_times_simple(double time_start, double delta, double time_end,
                            std::vector<double> &result) const;

  virtual bool process_ys(double &result);
  virtual bool process_y(double &result);
  virtual bool process_ye(double &result);

  virtual void compute_times(double time_start, double delta, double time_end,
                             const std::string &keyword,
                             std::vector<double> &result) const;

  virtual void parse_date(const std::string &spec, double *result) const;

  virtual void parse_interval_length(const std::string &spec, std::string &keyword,
                                     double *result) const;

protected:
  const Config::ConstPtr m_config;
  const units::System::Ptr m_unit_system;
  units::Unit m_time_units;
  double m_year_length;      //!< number of seconds in a year, for "mod" and "year fraction"
  double m_time_in_seconds, //!< current time, in seconds since the reference time
    m_run_start,                  //!< run start time, in seconds since the reference time
    m_run_end;                    //!< run end tim, in seconds since the reference time
  std::string m_calendar_string;       //!< CF calendar string
};

std::string reference_date_from_file(const PIO &nc,
                                     const std::string &time_name);

//! Create a Time instance by processing command-line options.
Time::Ptr time_from_options(MPI_Comm com, Config::ConstPtr config, units::System::Ptr system);

//! Initialize time from command-line options or from and input file (set using the `-i` option).
void initialize_time(MPI_Comm com, const std::string &dimension_name,
                     const Logger &log, Time &time);

} // end of namespace pism

#endif /* _PISMTIME_H_ */
