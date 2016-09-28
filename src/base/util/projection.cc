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
#include <cmath>                // fabs

#include "projection.hh"
#include "VariableMetadata.hh"
#include "error_handling.hh"
#include "io/PIO.hh"
#include "io/io_helpers.hh"

namespace pism {

MappingInfo::MappingInfo(const std::string &mapping_name, units::System::Ptr unit_system)
  : mapping(mapping_name, unit_system) {
  // empty
}

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

void check_consistency_epsg(const MappingInfo &info) {

  VariableMetadata epsg_mapping = epsg_to_cf(info.mapping.unit_system(), info.proj4);

  bool mapping_is_empty      = not info.mapping.has_attributes();
  bool epsg_mapping_is_empty = not epsg_mapping.has_attributes();

  if (mapping_is_empty and epsg_mapping_is_empty) {
    // empty mapping variables are equivalent
    return;
  } else {
    // Check if the "info.mapping" variable in the input file matches the EPSG code.
    // Check strings.
    VariableMetadata::StringAttrs strings = epsg_mapping.get_all_strings();
    VariableMetadata::StringAttrs::const_iterator j;
    for (j = strings.begin(); j != strings.end(); ++j) {
      if (not info.mapping.has_attribute(j->first)) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "PROJ.4 string \"%s\" requires %s = \"%s\",\n"
                                      "but the mapping variable has no %s.",
                                      info.proj4.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str());
      }

      std::string string = info.mapping.get_string(j->first);

      if (not (string == j->second)) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has %s = \"%s\".",
                                      info.proj4.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str(),
                                      string.c_str());
      }
    }

    // Check doubles
    VariableMetadata::DoubleAttrs doubles = epsg_mapping.get_all_doubles();
    VariableMetadata::DoubleAttrs::const_iterator k;
    for (k = doubles.begin(); k != doubles.end(); ++k) {
      if (not info.mapping.has_attribute(k->first)) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has no %s.",
                                      info.proj4.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str());
      }

      double value = info.mapping.get_double(k->first);

      if (fabs(value - k->second[0]) > 1e-12) {
        throw RuntimeError::formatted("inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has %s = %f.",
                                      info.proj4.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str(),
                                      value);
      }
    }
  }
}

MappingInfo get_projection_info(const PIO &input_file, const std::string &mapping_name,
                                units::System::Ptr unit_system) {
  MappingInfo result(mapping_name, unit_system);

  result.proj4 = input_file.get_att_text("PISM_GLOBAL", "proj4");

  std::string::size_type position = result.proj4.find("+init=epsg:");
  bool proj4_is_epsg = position != std::string::npos;

  if (input_file.inq_var(mapping_name)) {
    // input file has a mapping variable

    io::read_attributes(input_file, mapping_name, result.mapping);

    if (proj4_is_epsg) {
      // check consistency
      try {
        check_consistency_epsg(result);
      } catch (RuntimeError &e) {
        e.add_context("getting projection info from %s", input_file.inq_filename().c_str());
        throw;
      }
    } else {
      // use mapping read from input_file (can't check consistency here)
    }
  } else {
    // no mapping variable in the input file

    if (proj4_is_epsg) {
      result.mapping = epsg_to_cf(unit_system, result.proj4);
    } else {
      // leave mapping empty
    }
  }
  return result;
}

} // end of namespace pism
