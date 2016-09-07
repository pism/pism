/* Copyright (C) 2016 PISM Authors
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

#include "FrontalMelt.hh"

namespace pism {

FrontalMelt::FrontalMelt(IceGrid::ConstPtr g)
  : CalvingFrontRetreat(g, 1) {
  // empty
}

FrontalMelt::~FrontalMelt() {

}

void FrontalMelt::init() {
  m_log->message(2,
                 "* Initializing the parameterization of frontal melt...\n");

  // FIXME
}

void FrontalMelt::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                       std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
  (void) dict;
  (void) ts_dict;
}

void FrontalMelt::write_variables_impl(const std::set<std::string> &vars, const PIO& nc) {
  (void) vars;
  (void) nc;
}

void FrontalMelt::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  (void) keyword;
  (void) result;
}

void FrontalMelt::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                        IO_Type nctype) {
  (void) vars;
  (void) nc;
  (void) nctype;
}

void FrontalMelt::compute_calving_rate(const IceModelVec2CellType &mask,
                                       IceModelVec2S &result) {
  (void) mask;
  (void) result;
}

} // end of namespace pism
