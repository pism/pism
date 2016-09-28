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

#include <cstdlib>              // strtol

#include "projection.hh"
#include "VariableMetadata.hh"
#include "error_handling.hh"

namespace pism {

//! @brief Return CF-Convention "mapping" variable corresponding to an EPSG code specified in a
//! PROJ.4 string.
VariableMetadata epsg_to_cf(units::System::Ptr system, const std::string &proj4_string) {
  VariableMetadata mapping("mapping", system);

  std::string option = "+init=epsg:";
  std::string::size_type pos = proj4_string.find(option);
  std::string epsg_code = proj4_string.substr(pos + option.size());
  char *endptr = NULL;
  const char *str = epsg_code.c_str();
  long int epsg = strtol(str, &endptr, 10);
  if (endptr == str) {
    throw RuntimeError::formatted("failed to parse EPSG code at '%s' in PROJ.4 string '%s'",
                                  epsg_code.c_str(), proj4_string.c_str());
  }

  switch (epsg) {
  case 3413:
    mapping.set_double("latitude_of_projection_origin", 90.0);
    mapping.set_double("scale_factor_at_projection_origin", 1.0);
    mapping.set_double("straight_vertical_longitude_from_pole", -45.0);
    mapping.set_double("standard_parallel", 70.0);
    mapping.set_double("false_northing", 0.0);
    mapping.set_string("grid_mapping_name", "polar_stereographic");
    mapping.set_double("false_easting", 0.0);    
    break;
  case 3031:
    mapping.set_double("latitude_of_projection_origin", -90.0);
    mapping.set_double("scale_factor_at_projection_origin", 1.0);
    mapping.set_double("straight_vertical_longitude_from_pole", 0.0);
    mapping.set_double("standard_parallel", -71.0);
    mapping.set_double("false_northing", 0.0);
    mapping.set_string("grid_mapping_name", "polar_stereographic");
    mapping.set_double("false_easting", 0.0);
    break;
  default:
    throw RuntimeError::formatted("unknown EPSG code '%d' in PROJ.4 string '%s'",
                                  (int)epsg, proj4_string.c_str());
  }
  return mapping;
}

void check_mapping_equivalence(const VariableMetadata &mapping,
                               const std::string &proj4_string) {

  VariableMetadata epsg_mapping = epsg_to_cf(mapping.unit_system(), proj4_string);

  bool mapping_is_empty = (mapping.get_all_strings().empty() and
                           mapping.get_all_doubles().empty());

  bool epsg_mapping_is_empty = (epsg_mapping.get_all_strings().empty() and
                                epsg_mapping.get_all_doubles().empty());

  if (mapping_is_empty and epsg_mapping_is_empty) {
    // empty mapping variables are equivalent
    return;
  } else {
    // Check if the "mapping" variable in the input file matches the EPSG code.
    // Check strings.
    VariableMetadata::StringAttrs strings = epsg_mapping.get_all_strings();
    VariableMetadata::StringAttrs::const_iterator j;
    for (j = strings.begin(); j != strings.end(); ++j) {
      if (not mapping.has_attribute(j->first)) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has no %s.",
                                      proj4_string.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str());
      }

      if (not (mapping.get_string(j->first) == j->second)) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has %s = \"%s\".",
                                      proj4_string.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str(),
                                      mapping.get_string(j->first).c_str());
      }
    }

    // Check doubles
    VariableMetadata::DoubleAttrs doubles = epsg_mapping.get_all_doubles();
    VariableMetadata::DoubleAttrs::const_iterator k;
    for (k = doubles.begin(); k != doubles.end(); ++k) {
      if (not mapping.has_attribute(k->first)) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has no %s.",
                                      proj4_string.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str());
      }

      if (fabs(mapping.get_double(k->first) - k->second[0]) > 1e-12) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has %s = %f.",
                                      proj4_string.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str(),
                                      mapping.get_double(k->first));
      }
    }
  }
}

} // end of namespace pism
