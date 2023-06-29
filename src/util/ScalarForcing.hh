// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023 PISM Authors
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

#ifndef PISM_SCALARFORCING_H
#define PISM_SCALARFORCING_H

#include <memory>               // std::unique_ptr
#include <vector>               // std::vector

namespace pism {

class Context;
class File;
class Logger;

/*!
 * This class helps with loading and using scalar forcings such as scalar temperature
 * offsets.
 *
 * It processes command-line options, reads data from a file, and gets data corresponding
 * to a time interval [t, t+dt].
 */
class ScalarForcing {
public:
  ScalarForcing(const Context &ctx,
                const std::string &option_prefix,
                const std::string &variable_name,
                const std::string &units,
                const std::string &glaciological_units,
                const std::string &long_name);

  ScalarForcing(const Context &ctx,
                const std::string &filename,
                const std::string &variable_name,
                const std::string &units,
                const std::string &glaciological_units,
                const std::string &long_name,
                bool periodic);

  ~ScalarForcing();

  double value(double t) const;

  double average(double t, double dt) const;

private:
  // disable copy constructor and the assignment operator:
  ScalarForcing(const ScalarForcing &other);
  ScalarForcing& operator=(const ScalarForcing&);

  void initialize(const Context &ctx,
                  const std::string &filename,
                  const std::string &variable_name,
                  const std::string &units,
                  const std::string &glaciological_units,
                  const std::string &long_name,
                  bool periodic);

  double integral(double a, double b) const;

  struct Impl;
  Impl *m_impl;
};

} // end of namespace pism

#endif /* PISM_SCALARFORCING_H */
