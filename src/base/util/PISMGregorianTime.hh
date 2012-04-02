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

#include "PISMTime.hh"

class PISMGregorianTime : public PISMTime
{
public:
  PISMGregorianTime(MPI_Comm c, const NCConfigVariable &conf);
  virtual ~PISMGregorianTime() {}

  virtual PetscErrorCode init();

  virtual PetscErrorCode init_from_file(string filename);

  virtual double mod(double time, double period);

  virtual double year_fraction(double T);

  virtual string date(double T);

  virtual string date();

  virtual string start_date();

  virtual string end_date();

  virtual string units()
  { return CF_units(); }

  virtual bool use_reference_date()
  { return true; }

protected:
  utUnit ut_units;
};


#endif /* _PISMGREGORIANTIME_H_ */
