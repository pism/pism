// Copyright (C) 2012 PISM Authors
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

#ifndef _PISMGREGORIANTIME_H_
#define _PISMGREGORIANTIME_H_

class PISMGregorianTime
{
public:
  PISMGregorianTime(MPI_Comm c, const NCConfigVariable &conf);
  virtual ~PISMGregorianTime() {}

  //! \brief Intialize using command-line options.
  virtual PetscErrorCode init();

  //! \brief Returns time since the origin modulo period.
  virtual PetscReal mod(PetscReal time, PetscReal period)
  {
    if (period <= 0)
      return time;

    PetscReal tmp = time - floor(time / period) * period;

    if (fabs(tmp - period) < 1)
      tmp = 0;

    return tmp;
  }

  //! \brief Returns the fraction of a year passed since the last beginning of
  //! a year.
  virtual PetscReal year_fraction(PetscReal T)
  { return seconds_to_years(T) - floor(seconds_to_years(T)); }

  //! \brief Returns the year corresponding to time T. Only for reporting.
  virtual PetscReal year(PetscReal T)
  { return seconds_to_years(T); }

  //! \brief Returns current time, in years. Only for reporting.
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
  PetscErrorCode interval_to_seconds(string interval);
};


#endif /* _PISMGREGORIANTIME_H_ */
