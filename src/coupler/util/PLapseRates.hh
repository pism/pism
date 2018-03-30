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

#ifndef _PLAPSERATES_H_
#define _PLAPSERATES_H_

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"

namespace pism {

class Geometry;

template <class Model>
class PLapseRates : public Model {
public:
  PLapseRates(IceGrid::ConstPtr g, std::shared_ptr<Model> in)
    : Model(g, in) {
    m_temp_lapse_rate = 0.0;
  }

  virtual ~PLapseRates() {
    // empty
  }

protected:

protected:
  unsigned int m_bc_period;
  double m_bc_reference_time,          // in seconds
    m_temp_lapse_rate;
  std::string m_option_prefix;
};


} // end of namespace pism

#endif /* _PLAPSERATES_H_ */
