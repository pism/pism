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

#ifndef _POANOMALY_H_
#define _POANOMALY_H_

#include "pism/coupler/OceanModel.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {
namespace ocean {

//! @brief Reads and uses shelf_basal_mass_flux *anomalies* from a file
class Anomaly : public OceanModel
{
public:
  Anomaly(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in);
  virtual ~Anomaly();
protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual const IceModelVec2S& shelf_base_mass_flux_impl() const;

protected:
  IceModelVec2S::Ptr m_shelf_base_mass_flux;

  IceModelVec2T::Ptr m_shelf_base_mass_flux_anomaly;

};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POANOMALY_H_ */
