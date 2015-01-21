/* Copyright (C) 2015 PISM Authors
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

namespace pism {

Diagnostic::Diagnostic(const IceGrid &g)
  : m_grid(g) {
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
NCSpatialVariable Diagnostic::get_metadata(int N) {
  if (N >= m_dof) {
    return NCSpatialVariable(m_grid.config.get_unit_system(), "missing", m_grid);
  }

  return m_vars[N];
}

//! Define NetCDF variables corresponding to a diagnostic quantity.
void Diagnostic::define(const PIO &nc) {
  for (int j = 0; j < m_dof; ++j) {
    m_vars[j].define(nc, m_output_datatype, true);
  }
}

//! \brief A method for setting common variable attributes.
void Diagnostic::set_attrs(const std::string &my_long_name,
                           const std::string &my_standard_name,
                           const std::string &my_units,
                           const std::string &my_glaciological_units,
                           int N) {
  if (N >= m_dof) {
    throw RuntimeError::formatted("N (%d) >= m_dof (%d)", N, m_dof);
  }

  m_vars[N].set_string("pism_intent", "diagnostic");

  m_vars[N].set_string("long_name", my_long_name);

  m_vars[N].set_string("standard_name", my_standard_name);

  m_vars[N].set_units(my_units);

  if (not my_glaciological_units.empty()) {
    m_vars[N].set_glaciological_units(my_glaciological_units);
  }
}


} // end of namespace pism
