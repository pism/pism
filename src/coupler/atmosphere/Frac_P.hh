// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PAFPFORCING_H_
#define _PAFPFORCING_H_

#include "pism/coupler/AtmosphereModel.hh"

namespace pism {

class ScalarForcing;

namespace atmosphere {

class Frac_P : public AtmosphereModel {
public:
  Frac_P(IceGrid::ConstPtr g, std::shared_ptr<AtmosphereModel> in);
  virtual ~Frac_P();

private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  void init_timeseries_impl(const std::vector<double> &ts) const;

  const IceModelVec2S& mean_precipitation_impl() const;

  void precip_time_series_impl(int i, int j, std::vector<double> &values) const;

  mutable std::vector<double> m_offset_values;

  std::unique_ptr<ScalarForcing> m_forcing;

  IceModelVec2S::Ptr m_precipitation;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAFPFORCING_H_ */
