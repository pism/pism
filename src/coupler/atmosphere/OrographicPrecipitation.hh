// Copyright (C) 2018 PISM Authors
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
#include "pism/util/iceModelVec2T.hh"

namespace pism {

class Geometry;

namespace atmosphere {

class OrographicPrecipitationSerial;

class OrographicPrecipitation : public AtmosphereModel
{
public:
  OrographicPrecipitation(IceGrid::ConstPtr g);
  virtual ~OrographicPrecipitation();
private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& mean_annual_temp_impl() const;

  void begin_pointwise_access_impl() const;
  void end_pointwise_access_impl() const;

  void init_timeseries_impl(const std::vector<double> &ts) const;
  void temp_time_series_impl(int i, int j, std::vector<double> &values) const;

  IceModelVec2T::Ptr m_air_temp;
protected:
  std::string m_reference;

  IceModelVec2S m_precipitation;

  //! Storage on rank zero. Used to pass the load to the serial deformation model and get
  //! bed displacement back.
  petsc::Vec::Ptr m_work0;

  //! Ice-equivalent load thickness.
  IceModelVec2S m_surface;

  //! extended grid for the LT Model
  IceGrid::Ptr m_extended_grid;
  
  //! Serial orographic precipitation model.
  std::unique_ptr<OrographicPrecipitationSerial> m_serial_model;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAOROGRAPHICPRECIPITATION_H_ */
