// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <memory>               // std::unique_ptr

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/Timeseries.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"

namespace pism {

template<class Model, class Mod>
class PScalarForcing : public Mod
{
public:
  PScalarForcing(IceGrid::ConstPtr g, std::shared_ptr<Model> in)
    : Mod(g, in), m_current_forcing(NAN) {
    // empty
  }

  virtual ~PScalarForcing()
  {
    // empty
  }

protected:
  virtual void update_impl(double t, double dt)
  {
    Mod::m_t  = Mod::m_grid->ctx()->time()->mod(t - m_bc_reference_time, m_bc_period);
    Mod::m_dt = dt;

    // this call will update the input model, too
    Mod::update_impl(t, dt);

    m_current_forcing = (*m_offset)(Mod::m_t + 0.5*Mod::m_dt);
  }

  virtual void init_internal()
  {
    IceGrid::ConstPtr g = Mod::m_grid;

    options::String file(m_option_prefix + "_file", "Specifies a file with scalar offsets");
    options::Integer period(m_option_prefix + "_period",
                            "Specifies the length of the climate data period", 0);
    options::Real bc_reference_year(m_option_prefix + "_reference_year",
                                    "Boundary condition reference year", 0.0);

    if (not file.is_set()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "command-line option %s_file is required.",
                                    m_option_prefix.c_str());
    }

    if (period.value() < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid %s_period %d (period length cannot be negative)",
                                    m_option_prefix.c_str(), period.value());
    }
    m_bc_period = (unsigned int)period;

    if (bc_reference_year.is_set()) {
      m_bc_reference_time = units::convert(Mod::m_sys, bc_reference_year, "years", "seconds");
    } else {
      m_bc_reference_time = 0;
    }

    Mod::m_log->message(2,
                        "  reading %s data from forcing file %s...\n",
                        m_offset->name().c_str(), file->c_str());

    PIO nc(g->com, "netcdf3", file, PISM_READONLY);
    {
      m_offset->read(nc, *g->ctx()->time(), *g->ctx()->log());
    }
  }

  //! Apply the current forcing as an offset.
  void offset_data(IceModelVec2S &result) const {
    result.shift(m_current_forcing);
  }

  //! Apply the current forcing as a scaling factor.
  void scale_data(IceModelVec2S &result) const {
    result.scale(m_current_forcing);
  }

  std::unique_ptr<Timeseries> m_offset;
  std::string m_filename, m_offset_name, m_option_prefix;

  unsigned int m_bc_period;       // in years
  double m_bc_reference_time;  // in seconds
  // current scalar forcing
  double m_current_forcing;
};


} // end of namespace pism

#endif /* _PSCALARFORCING_H_ */
