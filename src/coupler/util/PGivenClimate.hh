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

#ifndef _PGIVENCLIMATE_H_
#define _PGIVENCLIMATE_H_

#include <memory>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Component.hh"

namespace pism {

class Geometry;

template <class Model>
class PGivenClimate : public Model
{
public:
  PGivenClimate(IceGrid::ConstPtr g, std::shared_ptr<Model> in)
    : Model(g, in) {}

  virtual ~PGivenClimate() {
    // empty
  }

protected:


protected:
  std::string m_filename;

  unsigned int m_bc_period;       // in (integer) years
  double m_bc_reference_time;  // in seconds
};

} // end of namespace pism

#endif /* _PGIVENCLIMATE_H_ */
