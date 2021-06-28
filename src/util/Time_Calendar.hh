// Copyright (C) 2012, 2013, 2014, 2015, 2017, 2019, 2021 PISM Authors
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

#include "Time.hh"
#include "Units.hh"

namespace pism {

class Time_Calendar : public Time
{
public:
  Time_Calendar(MPI_Comm c, Config::ConstPtr conf,
                const std::string &calendar,
                units::System::Ptr units_system);
  virtual ~Time_Calendar() = default;

  virtual void init(const Logger &log);

  virtual void init_from_file(const std::string &filename, const Logger &log,
                              bool set_start_time);

  virtual void init_from_input_file(const File &nc,
                                    const std::string &time_name,
                                    const Logger &log);

  virtual double mod(double time, unsigned int) const;

  virtual double year_fraction(double T) const;

  using Time::date;
  virtual std::string date(double T) const;

  virtual std::string units_string() const {
    return CF_units_string();
  }

  virtual std::string CF_units_string() const {
    return m_time_units.format();
  }

  virtual std::string CF_units_to_PISM_units(const std::string &input) const {
    return input;               // return unchanged CF units
  }

  virtual bool use_reference_date() {
    return true;
  }

  virtual double calendar_year_start(double T) const;

  virtual double increment_date(double T, int years) const;

protected:
  virtual void compute_times(double time_start, double time_end,
                             const Interval &interval,
                             std::vector<double> &result) const;

  virtual bool process_ys(double &result);
  virtual bool process_y(double &result);
  virtual bool process_ye(double &result);

  virtual double parse_date(const std::string &spec) const;

  virtual Interval parse_interval_length(const std::string &spec) const;

  void compute_times_monthly(std::vector<double> &result) const;

  void compute_times_yearly(std::vector<double> &result) const;
private:
  MPI_Comm m_com;
  bool m_simple_calendar;
  // Hide copy constructor / assignment operator.
  Time_Calendar(Time_Calendar const &);
  Time_Calendar & operator=(Time_Calendar const &);
};


} // end of namespace pism

#endif /* _PISMGREGORIANTIME_H_ */
