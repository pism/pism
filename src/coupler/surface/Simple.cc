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

#include <cassert>
#include <gsl/gsl_math.h>

#include "Simple.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace surface {

///// Simple PISM surface model.
Simple::Simple(IceGrid::ConstPtr g)
  : SurfaceModel(g) {
}

void Simple::init_impl(const Geometry &geometry) {
  assert(m_atmosphere != NULL);
  m_atmosphere->init(geometry);

  m_log->message(2,
             "* Initializing the simplest PISM surface (snow) processes model Simple.\n"
             "  It passes atmospheric state directly to upper ice fluid surface:\n"
             "    surface mass balance          := precipitation,\n"
             "    ice upper surface temperature := 2m air temperature.\n");
}

MaxTimestep Simple::max_timestep_impl(double t) const {
  (void) t;
  if (m_atmosphere) {
    return m_atmosphere->max_timestep(t);
  } else {
    return MaxTimestep("surface simple");
  }
}

void Simple::update_impl(const Geometry &geometry, double t, double dt) {
  m_t = t;
  m_dt = dt;
  if (m_atmosphere) {
    m_atmosphere->update(geometry, t, dt);
  }
}

void Simple::mass_flux_impl(IceModelVec2S &result) const {
  m_atmosphere->mean_precipitation(result);
}

void Simple::temperature_impl(IceModelVec2S &result) const {
  m_atmosphere->mean_annual_temp(result);
}

} // end of namespace surface
} // end of namespace pism
