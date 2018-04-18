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

#ifndef _SCALARFORCING_H_
#define _SCALARFORCING_H_

#include <memory>               // std::unique_ptr

#include "pism/util/Context.hh"

namespace pism {

class Timeseries;

/*!
 * This class helps with loading and using scalar forcings such as scalar temperature
 * offsets.
 *
 * It processes command-line options, reads data from a file, and gets data corresponding
 * to a time interval [t, t+dt].
 */
class ScalarForcing {
public:
  ScalarForcing(Context::ConstPtr ctx,
                const std::string &option_prefix,
                const std::string &offset_name,
                const std::string &units,
                const std::string &glaciological_units,
                const std::string &long_name);
  ~ScalarForcing();

  void init();
  void update(double t, double dt);

  double value() const;
  double value(double t) const;
protected:
  Context::ConstPtr m_ctx;

  std::unique_ptr<Timeseries> m_data;

  std::string m_option_prefix;

  // in years
  unsigned int m_bc_period;
  // in seconds
  double m_bc_reference_time;

  double m_current;
};


} // end of namespace pism

#endif /* _SCALARFORCING_H_ */
