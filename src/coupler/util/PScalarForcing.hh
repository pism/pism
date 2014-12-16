// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PSCALARFORCING_H_
#define _PSCALARFORCING_H_

#include "IceGrid.hh"
#include "iceModelVec.hh"
#include "Timeseries.hh"
#include "pism_options.hh"
#include "PISMTime.hh"

#include "error_handling.hh"
#include "PISMConfig.hh"

namespace pism {

template<class Model, class Mod>
class PScalarForcing : public Mod
{
public:
  PScalarForcing(IceGrid &g, Model* in)
    : Mod(g, in), input(in) {}

  virtual ~PScalarForcing()
  {
    if (offset) {
      delete offset;
    }
  }

  virtual void update(double my_t, double my_dt)
  {
    Mod::m_t  = Mod::grid.time->mod(my_t - bc_reference_time, bc_period);
    Mod::m_dt = my_dt;

    Mod::input_model->update(my_t, my_dt);
  }

protected:
  virtual void init_internal()
  {
    bool file_set, bc_period_set, bc_ref_year_set;

    IceGrid &g = Mod::grid;

    double bc_period_years = 0,
      bc_reference_year = 0;

    {
      OptionsString(option_prefix + "_file", "Specifies a file with scalar offsets",
                    filename, file_set);
      OptionsReal(option_prefix + "_period", "Specifies the length of the climate data period",
                  bc_period_years, bc_period_set);
      OptionsReal(option_prefix + "_reference_year", "Boundary condition reference year",
                  bc_reference_year, bc_ref_year_set);
    }

    if (file_set == false) {
      throw RuntimeError::formatted("command-line option %s_file is required.",
                                    option_prefix.c_str());
    }

    if (bc_period_set) {
      bc_period = (unsigned int)bc_period_years;
    } else {
      bc_period = 0;
    }

    if (bc_ref_year_set) {
      bc_reference_time = g.convert(bc_reference_year, "years", "seconds");
    } else {
      bc_reference_time = 0;
    }


    verbPrintf(2, g.com,
               "  reading %s data from forcing file %s...\n",
               offset->short_name.c_str(), filename.c_str());

    PIO nc(g.com, "netcdf3", g.config.get_unit_system());
    nc.open(filename, PISM_READONLY);
    {
      offset->read(nc, g.time);
    }
    nc.close();
  }

  //! Apply offset as an offset
  void offset_data(IceModelVec2S &result) {
    result.shift((*offset)(Mod::m_t + 0.5*Mod::m_dt));
  }

  //! Apply offset as a scaling factor
  void scale_data(IceModelVec2S &result) {
    result.scale((*offset)(Mod::m_t + 0.5*Mod::m_dt));
  }

  Model *input;
  Timeseries *offset;
  std::string filename, offset_name, option_prefix;

  unsigned int bc_period;       // in years
  double bc_reference_time;  // in seconds
};


} // end of namespace pism

#endif /* _PSCALARFORCING_H_ */
