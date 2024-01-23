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

#ifndef _CLIMATEINDEXWEIGHTS_H_
#define _CLIMATEINDEXWEIGHTS_H_

#include "pism/util/ScalarForcing.hh"
#include "pism/util/Context.hh"

namespace pism {

class ScalarForcing;

class ClimateIndexWeights {
public:
    ClimateIndexWeights(const Context &ctx);
    virtual ~ClimateIndexWeights() = default;

   // void init_weights();
    void update_weights(double t, double dt, double &m_W0, double &m_W1, double &m_W1X);
    // const std::vector<double> get_weights() const;

protected:
    // glacial index timeseries
    std::unique_ptr<ScalarForcing> m_index;

    // glacial index weights
    double m_current, m_ref_value, m_max_value;//, m_W0, m_W1, m_W1X;
};

} // end of namespace pism

#endif /* _CLIMATEINDEXWEIGHTS_H_ */