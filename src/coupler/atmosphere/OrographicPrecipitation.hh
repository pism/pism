// Copyright (C) 2018, 2020, 2021 PISM Authors
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

#ifndef _PAOROGRAPHICPRECIPITATION_H_
#define _PAOROGRAPHICPRECIPITATION_H_

#include "pism/coupler/AtmosphereModel.hh"

namespace pism {

class Geometry;

namespace atmosphere {

class OrographicPrecipitationSerial;

class OrographicPrecipitation : public AtmosphereModel {
public:
  OrographicPrecipitation(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~OrographicPrecipitation();

private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const IceModelVec2S &mean_precipitation_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;

protected:
  std::string m_reference;

  IceModelVec2S::Ptr m_precipitation;

  //! Storage on rank zero. Used to pass the load to the serial orographic precipitation
  //! Model.
  std::shared_ptr<petsc::Vec> m_work0;

  //! Serial orographic precipitation model.
  std::unique_ptr<OrographicPrecipitationSerial> m_serial_model;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAOROGRAPHICPRECIPITATION_H_ */
