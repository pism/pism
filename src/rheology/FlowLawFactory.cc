// Copyright (C) 2009--2018, 2023, 2025 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cassert>
#include <string>

#include "pism/rheology/FlowLawFactory.hh"
#include "pism/util/Config.hh"
#include "pism/util/error_handling.hh"

#include "pism/rheology/IsothermalGlen.hh"
#include "pism/rheology/PatersonBudd.hh"
#include "pism/rheology/GPBLD.hh"
#include "pism/rheology/Hooke.hh"
#include "pism/rheology/PatersonBuddCold.hh"
#include "pism/rheology/PatersonBuddWarm.hh"
#include "pism/rheology/GoldsbyKohlstedt.hh"

namespace pism {
namespace rheology {

using ECPtr = std::shared_ptr<EnthalpyConverter>;

FlowLaw *create_isothermal_glen(double exponent, const Config &config, ECPtr EC) {
  return new (IsothermalGlen)(exponent, config, EC);
}

FlowLaw *create_pb(double exponent, const Config &config, ECPtr EC) {
  return new (PatersonBudd)(exponent, config, EC);
}

FlowLaw *create_gpbld(double exponent, const Config &config, ECPtr EC) {
  return new (GPBLD)(exponent, config, EC);
}

FlowLaw *create_hooke(double exponent, const Config &config, ECPtr EC) {
  return new (Hooke)(exponent, config, EC);
}

FlowLaw *create_arr(double exponent, const Config &config, ECPtr EC) {
  return new (PatersonBuddCold)(exponent, config, EC);
}

FlowLaw *create_arrwarm(double exponent, const Config &config, ECPtr EC) {
  return new (PatersonBuddWarm)(exponent, config, EC);
}

FlowLaw *create_goldsby_kohlstedt(double exponent, const Config &config, ECPtr EC) {
  return new (GoldsbyKohlstedt)(exponent, config, EC);
}

FlowLawFactory::FlowLawFactory(std::shared_ptr<const Config> conf,
                               std::shared_ptr<EnthalpyConverter> EC)
  : m_config(conf), m_EC(EC) {

  m_flow_laws = { { ICE_ISOTHERMAL_GLEN, create_isothermal_glen },
                  { ICE_PB, create_pb },
                  { ICE_GPBLD, create_gpbld },
                  { ICE_HOOKE, create_hooke },
                  { ICE_ARR, create_arr },
                  { ICE_ARRWARM, create_arrwarm },
                  { ICE_GOLDSBY_KOHLSTEDT, create_goldsby_kohlstedt } };
}

void FlowLawFactory::add(const std::string &name, FlowLawCreator icreate) {
  m_flow_laws[name] = icreate;
}

void FlowLawFactory::remove(const std::string &name) {
  m_flow_laws.erase(name);
}

std::shared_ptr<FlowLaw> FlowLawFactory::create(const std::string &type_name, double exponent) {
  // find the function that can create selected flow law:
  if (m_flow_laws[type_name] == nullptr) {
    throw RuntimeError::formatted(
        PISM_ERROR_LOCATION, "Selected ice flow law \"%s\" is not available", type_name.c_str());
  }

  auto r = m_flow_laws[type_name];

  // create an FlowLaw instance:
  return std::shared_ptr<FlowLaw>((*r)(exponent, *m_config, m_EC));
}

} // end of namespace rheology
} // end of namespace pism
