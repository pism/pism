// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _PSCALARFORCING_H_
#define _PSCALARFORCING_H_

#include "IceGrid.hh"
#include "iceModelVec.hh"
#include "Timeseries.hh"
#include "pism_options.hh"
#include "PISMTime.hh"

template<class Model, class Mod>
class PScalarForcing : public Mod
{
public:
  PScalarForcing(IceGrid &g, const NCConfigVariable &conf, Model* in)
    : Mod(g, conf, in), input(in) {}
  virtual ~PScalarForcing()
  {
    if (offset)
      delete offset;
  }

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  {
    Mod::t  = Mod::grid.time->mod(my_t - bc_reference_time, bc_period);
    Mod::dt = my_dt;

    PetscErrorCode ierr = Mod::input_model->update(my_t, my_dt); CHKERRQ(ierr);
    return 0;
  }

protected:
  virtual PetscErrorCode init_internal()
  {
    PetscErrorCode ierr;
    bool file_set, bc_period_set, bc_ref_year_set;

    IceGrid &g = Mod::grid;

    PetscReal bc_period_years = 0,
      bc_reference_year = 0;

    ierr = PetscOptionsBegin(g.com, "", "Scalar forcing options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString(option_prefix + "_file", "Specifies a file with scalar offsets",
                               filename, file_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_period", "Specifies the length of the climate data period",
                             bc_period_years, bc_period_set); CHKERRQ(ierr);
      ierr = PISMOptionsReal(option_prefix + "_reference_year", "Boundary condition reference year",
                             bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (file_set == false) {
      ierr = verbPrintf(2, g.com, "  WARNING: %s_file is not set; forcing is inactive.\n",
                        option_prefix.c_str()); CHKERRQ(ierr);
      delete offset;
      offset = NULL;
    }

    if (bc_period_set) {
      bc_period = g.time->years_to_seconds(bc_period_years);
    } else {
      bc_period = 0;
    }

    if (bc_ref_year_set) {
      bc_reference_time = g.time->years_to_seconds(bc_reference_year);
    } else {
      bc_reference_time = 0;
    }

    if (offset) {
      ierr = verbPrintf(2, g.com,
                        "  reading %s data from forcing file %s...\n",
                        offset->short_name.c_str(), filename.c_str());
      CHKERRQ(ierr);

      ierr = offset->read(filename.c_str(), g.time->use_reference_date()); CHKERRQ(ierr);
    }

    return 0;
  }

  PetscErrorCode offset_data(IceModelVec2S &result) {
    if (offset) {
      PetscErrorCode ierr = result.shift((*offset)(Mod::t + 0.5*Mod::dt)); CHKERRQ(ierr);
    }
    return 0;
  }

  Model *input;
  Timeseries *offset;
  string filename, offset_name, option_prefix;

  PetscReal bc_period,          // in seconds
    bc_reference_time;          // in seconds

};


#endif /* _PSCALARFORCING_H_ */
