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

    PetscReal tmp = time - floor(time / period) * period;

    if (fabs(tmp - period) < 1)
      tmp = 0;

    return tmp;
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

  virtual PetscReal start_year()
  { return run_start / secpera; }

  virtual PetscReal end_year()
  { return run_end / secpera; }

};


#endif /* _PISMGREGORIANTIME_H_ */
