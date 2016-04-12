/* Copyright (C) 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMDiagnostic.hh"
#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "base/util/Logger.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

Diagnostic::Diagnostic(IceGrid::ConstPtr g)
  : m_grid(g), m_sys(g->ctx()->unit_system()), m_config(g->ctx()->config()) {
  m_output_datatype = PISM_FLOAT;
  m_dof = 1;
}

Diagnostic::~Diagnostic() {
  // empty
}

//! \brief Update a cumulative quantity needed to compute a rate of change.
//! So far we there is only one such quantity: the rate of change of the ice
//! thickness.
void Diagnostic::update_cumulative() {
  // the default implementation is empty
}

//! Get the number of NetCDF variables corresponding to a diagnostic quantity.
int Diagnostic::get_nvars() {
  return m_dof;
}

//! Reset vertical levels corresponding to the z dimension.
/** This is called after the automatic grid extension.
 */
void Diagnostic::set_zlevels(std::vector<double> &zlevels) {
  for (int j = 0; j < m_dof; ++j) {
    if (m_vars[j].get_z().get_name() == "z") {
      m_vars[j].set_levels(zlevels);
    }
  }
}

//! Get a metadata object corresponding to variable number N.
SpatialVariableMetadata Diagnostic::get_metadata(int N) {
  if (N >= m_dof) {
    return SpatialVariableMetadata(m_grid->ctx()->unit_system(), "missing");
  }

  return m_vars[N];
}

//! Define NetCDF variables corresponding to a diagnostic quantity.
void Diagnostic::define(const PIO &nc) {
  std::string order = m_grid->ctx()->config()->get_string("output_variable_order");
  for (int j = 0; j < m_dof; ++j) {
    io::define_spatial_variable(m_vars[j], *m_grid, nc, m_output_datatype, order, true);
  }
}

//! \brief A method for setting common variable attributes.
void Diagnostic::set_attrs(const std::string &long_name,
                           const std::string &standard_name,
                           const std::string &units,
                           const std::string &glaciological_units,
                           int N) {
  if (N >= m_dof) {
    throw RuntimeError::formatted("N (%d) >= m_dof (%d)", N, m_dof);
  }

  m_vars[N].set_string("pism_intent", "diagnostic");

  m_vars[N].set_string("long_name", long_name);

  m_vars[N].set_string("standard_name", standard_name);

  m_vars[N].set_string("units", units);

  if (not glaciological_units.empty()) {
    m_vars[N].set_string("glaciological_units", glaciological_units);
  }
}

IceModelVec::Ptr Diagnostic::compute() {
  // use the name of the first variable
  std::vector<std::string> names;
  for (unsigned int j = 0; j < m_vars.size(); ++j) {
    names.push_back(m_vars[j].get_name());
  }
  std::string all_names = join(names, ",");

  m_grid->ctx()->log()->message(3, "-  Computing %s...\n", all_names.c_str());
  IceModelVec::Ptr result = this->compute_impl();
  m_grid->ctx()->log()->message(3, "-  Done computing %s.\n", all_names.c_str());
  return result;
}

} // end of namespace pism
