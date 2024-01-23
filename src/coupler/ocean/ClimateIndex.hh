// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023 PISM Authors
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

#ifndef _CLIMATEINDEXOCEAN_H_
#define _CLIMATEINDEXOCEAN_H_

#include "CompleteOceanModel.hh"
#include "pism/coupler/util/ClimateIndexWeights.hh"

namespace pism {

class ClimateIndexWeights;

namespace ocean {

class ClimateIndex : public CompleteOceanModel
{
public:
    ClimateIndex(std::shared_ptr<const Grid> grid);
    virtual ~ClimateIndex() = default;
    virtual void init_forcing();
    virtual void update_forcing(double t, double dt, array::Scalar &theta_ocean, array::Scalar &salinity_ocean);

protected:
   std::unique_ptr<ClimateIndexWeights> m_climate_index;
   bool use_1X;
   double m_w0, m_w1, m_w1X;
   std::string m_reference;

  // ocean state at glacial index reference value (e.g. Present-Day, PD)
  array::Scalar m_theta_ocean_ref, m_salinity_ocean_ref;
  // ocean state anomaly at glacial index zero (Last Glacial Maximum, LGM)
  array::Scalar m_theta_ocean_anomaly_0, m_salinity_ocean_anomaly_0;
  // ocean state anomaly at glacial index one (e.g. the Eemian, LIG)
  array::Scalar m_theta_ocean_anomaly_1, m_salinity_ocean_anomaly_1;

  array::Scalar m_theta_ocean_anomaly_1X, m_salinity_ocean_anomaly_1X;

};

} // end of namespace ocean
} // end of namespace pism

#endif /* _CLIMATEINDEXOCEAN_H_ */
