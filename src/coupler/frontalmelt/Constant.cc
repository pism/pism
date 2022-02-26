// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2022 PISM Authors
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

#include "Constant.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace frontalmelt {

Constant::Constant(IceGrid::ConstPtr g)
  : FrontalMelt(g) {
  m_frontal_melt_rate = std::make_shared<array::Scalar>(g, "frontal_melt_rate");
  m_frontal_melt_rate->set_attrs("diagnostic", "frontal melt rate",
                                 "m s-1", "m day-1", "", 0);
}

void Constant::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {
  (void) t;
  (void) dt;

  const auto &cell_type = inputs.geometry->cell_type;

  const double
    melt_rate = m_config->get_number("frontal_melt.constant.melt_rate", "m second-1");

  IceModelVec::AccessList list{&cell_type, m_frontal_melt_rate.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (apply(cell_type, i, j)) {
      (*m_frontal_melt_rate)(i, j) = melt_rate;
    } else {
      (*m_frontal_melt_rate)(i, j) = 0.0;
    }
  }
}

const array::Scalar& Constant::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}
  
void Constant::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the constant frontal melt model...\n"
                 "  Frontal melt rate set to %f m/year.\n",
                 m_config->get_number("frontal_melt.constant.melt_rate", "m year-1"));
}

MaxTimestep Constant::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("frontal_melt constant");
}


} // end of namespape frontalmelt
} // end of namespace pism
