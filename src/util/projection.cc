/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020 PISM Authors
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
#include "io/File.hh"
#include "io/io_helpers.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"

#include "pism/pism_config.hh"

#if (Pism_USE_PROJ==1)
#include "pism/util/Proj.hh"
#endif

namespace pism {

MappingInfo::MappingInfo(const std::string &mapping_name, units::System::Ptr unit_system)
  : mapping(mapping_name, unit_system) {
  // empty
}

//! @brief Return CF-Convention "mapping" variable corresponding to an EPSG code specified in a
//! PROJ string.
VariableMetadata epsg_to_cf(units::System::Ptr system, const std::string &proj_string) {
  VariableMetadata mapping("mapping", system);

  int auth_len = 5;             // length of "epsg:"
  std::string::size_type position = std::string::npos;
  for (auto &auth : {"epsg:", "EPSG:"}) {
    position = proj_string.find(auth);
    if (position != std::string::npos) {
      break;
    }
  }

  if (position == std::string::npos) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "could not parse EPSG code '%s'", proj_string.c_str());
  }

  long int epsg = 0;
  {
    std::string epsg_code = proj_string.substr(position + auth_len);
    const char *str = epsg_code.c_str();
    char *endptr = NULL;
    epsg = strtol(str, &endptr, 10);
    if (endptr == str) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "failed to parse EPSG code at '%s' in PROJ string '%s'",
                                    epsg_code.c_str(), proj_string.c_str());
    }
  }

  switch (epsg) {
  case 3413:
    mapping.set_number("latitude_of_projection_origin", 90.0);
    mapping.set_number("scale_factor_at_projection_origin", 1.0);
    mapping.set_number("straight_vertical_longitude_from_pole", -45.0);
    mapping.set_number("standard_parallel", 70.0);
    mapping.set_number("false_northing", 0.0);
    mapping.set_string("grid_mapping_name", "polar_stereographic");
    mapping.set_number("false_easting", 0.0);
    break;
  case 3031:
    mapping.set_number("latitude_of_projection_origin", -90.0);
    mapping.set_number("scale_factor_at_projection_origin", 1.0);
    mapping.set_number("straight_vertical_longitude_from_pole", 0.0);
    mapping.set_number("standard_parallel", -71.0);
    mapping.set_number("false_northing", 0.0);
    mapping.set_string("grid_mapping_name", "polar_stereographic");
    mapping.set_number("false_easting", 0.0);
    break;
  case 3057:
    mapping.set_string("grid_mapping_name", "lambert_conformal_conic") ;
    mapping.set_number("longitude_of_central_meridian", -19.) ;
    mapping.set_number("false_easting", 500000.) ;
    mapping.set_number("false_northing", 500000.) ;
    mapping.set_number("latitude_of_projection_origin", 65.) ;
    mapping.set_numbers("standard_parallel", {64.25, 65.75}) ;
    mapping.set_string("long_name", "CRS definition") ;
    mapping.set_number("longitude_of_prime_meridian", 0.) ;
    mapping.set_number("semi_major_axis", 6378137.) ;
    mapping.set_number("inverse_flattening", 298.257222101) ;
    break;
  case 5936:
    mapping.set_number("latitude_of_projection_origin", 90.0);
    mapping.set_number("scale_factor_at_projection_origin", 1.0);
    mapping.set_number("straight_vertical_longitude_from_pole", -150.0);
    mapping.set_number("standard_parallel", 90.0);
    mapping.set_number("false_northing", 2000000.0);
    mapping.set_string("grid_mapping_name", "polar_stereographic");
    mapping.set_number("false_easting", 2000000.0);
    break;
  case 26710:
    mapping.set_number("longitude_of_central_meridian", -123.0);
    mapping.set_number("false_easting", 500000.0);
    mapping.set_number("false_northing", 0.0);
    mapping.set_string("grid_mapping_name", "transverse_mercator");
    mapping.set_number("inverse_flattening", 294.978698213898);
    mapping.set_number("latitude_of_projection_origin", 0.0);
    mapping.set_number("scale_factor_at_central_meridian", 0.9996);
    mapping.set_number("semi_major_axis", 6378206.4);
    mapping.set_string("unit", "metre");
    break;
  default:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "unknown EPSG code '%d' in PROJ string '%s'",
                                  (int)epsg, proj_string.c_str());
  }

  return mapping;
}

void check_consistency_epsg(const MappingInfo &info) {

  VariableMetadata epsg_mapping = epsg_to_cf(info.mapping.unit_system(), info.proj);

  bool mapping_is_empty      = not info.mapping.has_attributes();
  bool epsg_mapping_is_empty = not epsg_mapping.has_attributes();

  if (mapping_is_empty and epsg_mapping_is_empty) {
    // empty mapping variables are equivalent
    return;
  } else {
    // Check if the "info.mapping" variable in the input file matches the EPSG code.
    // Check strings.
    for (auto s : epsg_mapping.get_all_strings()) {
      if (not info.mapping.has_attribute(s.first)) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistent metadata:\n"
                                      "PROJ string \"%s\" requires %s = \"%s\",\n"
                                      "but the mapping variable has no %s.",
                                      info.proj.c_str(),
                                      s.first.c_str(), s.second.c_str(),
                                      s.first.c_str());
      }

      std::string string = info.mapping.get_string(s.first);

      if (not (string == s.second)) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has %s = \"%s\".",
                                      info.proj.c_str(),
                                      s.first.c_str(), s.second.c_str(),
                                      s.first.c_str(),
                                      string.c_str());
      }
    }

    // Check doubles
    for (auto d : epsg_mapping.get_all_doubles()) {
      if (not info.mapping.has_attribute(d.first)) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has no %s.",
                                      info.proj.c_str(),
                                      d.first.c_str(), d.second[0],
                                      d.first.c_str());
      }

      double value = info.mapping.get_number(d.first);

      if (fabs(value - d.second[0]) > 1e-12) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has %s = %f.",
                                      info.proj.c_str(),
                                      d.first.c_str(), d.second[0],
                                      d.first.c_str(),
                                      value);
      }
    }
  }
}

MappingInfo get_projection_info(const File &input_file, const std::string &mapping_name,
                                units::System::Ptr unit_system) {
  MappingInfo result(mapping_name, unit_system);

  result.proj = input_file.read_text_attribute("PISM_GLOBAL", "proj");

  bool proj_is_epsg = false;
  for (auto &auth : {"epsg:", "EPSG:"}) {
    if (result.proj.find(auth) != std::string::npos) {
      proj_is_epsg = true;
      break;
    }
  }

  if (input_file.find_variable(mapping_name)) {
    // input file has a mapping variable

    io::read_attributes(input_file, mapping_name, result.mapping);

    if (proj_is_epsg) {
      // check consistency
      try {
        check_consistency_epsg(result);
      } catch (RuntimeError &e) {
        e.add_context("getting projection info from %s", input_file.filename().c_str());
        throw;
      }
    } else {
      // use mapping read from input_file (can't check consistency here)
    }
  } else {
    // no mapping variable in the input file

    if (proj_is_epsg) {
      result.mapping = epsg_to_cf(unit_system, result.proj);
    } else {
      // leave mapping empty
    }
  }
  return result;
}

enum LonLat {LONGITUDE, LATITUDE};

#if (Pism_USE_PROJ==1)

//! Computes the area of a triangle using vector cross product.
static double triangle_area(double *A, double *B, double *C) {
  double V1[3], V2[3];
  for (int j = 0; j < 3; ++j) {
    V1[j] = B[j] - A[j];
    V2[j] = C[j] - A[j];
  }
  using std::pow;
  using std::sqrt;
  return 0.5*sqrt(pow(V1[1]*V2[2] - V2[1]*V1[2], 2) +
                  pow(V1[0]*V2[2] - V2[0]*V1[2], 2) +
                  pow(V1[0]*V2[1] - V2[0]*V1[1], 2));
}

void compute_cell_areas(const std::string &projection, IceModelVec2S &result) {
  IceGrid::ConstPtr grid = result.grid();

  Proj pism_to_geocent(projection, "+proj=geocent +datum=WGS84 +ellps=WGS84");

// Cell layout:
// (nw)        (ne)
// +-----------+
// |           |
// |           |
// |     o     |
// |           |
// |           |
// +-----------+
// (sw)        (se)

  double dx2 = 0.5 * grid->dx(), dy2 = 0.5 * grid->dy();

  IceModelVec::AccessList list(result);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      x = grid->x(i),
      y = grid->y(j);
    double
      x_nw = x - dx2, y_nw = y + dy2,
      x_ne = x + dx2, y_ne = y + dy2,
      x_se = x + dx2, y_se = y - dy2,
      x_sw = x - dx2, y_sw = y - dy2;

    PJ_COORD in, out;

    in.xy = {x_nw, y_nw};
    out = proj_trans(*pism_to_geocent, PJ_FWD, in);
    double nw[3] = {out.xyz.x, out.xyz.y, out.xyz.z};

    in.xy = {x_ne, y_ne};
    out = proj_trans(*pism_to_geocent, PJ_FWD, in);
    double ne[3] = {out.xyz.x, out.xyz.y, out.xyz.z};

    in.xy = {x_se, y_se};
    out = proj_trans(*pism_to_geocent, PJ_FWD, in);
    double se[3] = {out.xyz.x, out.xyz.y, out.xyz.z};

    in.xy = {x_sw, y_sw};
    out = proj_trans(*pism_to_geocent, PJ_FWD, in);
    double sw[3] = {out.xyz.x, out.xyz.y, out.xyz.z};

    result(i, j) = triangle_area(sw, se, ne) + triangle_area(ne, nw, sw);
  }
}

static void compute_lon_lat(const std::string &projection,
                            LonLat which, IceModelVec2S &result) {

  Proj crs(projection, "EPSG:4326");

// Cell layout:
// (nw)        (ne)
// +-----------+
// |           |
// |           |
// |     o     |
// |           |
// |           |
// +-----------+
// (sw)        (se)

  IceGrid::ConstPtr grid = result.grid();

  IceModelVec::AccessList list{&result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    PJ_COORD in, out;

    in.xy = {grid->x(i), grid->y(j)};
    out = proj_trans(*crs, PJ_FWD, in);

    if (which == LONGITUDE) {
      result(i, j) = out.lp.phi;
    } else {
      result(i, j) = out.lp.lam;
    }
  }
}

static void compute_lon_lat_bounds(const std::string &projection,
                                   LonLat which,
                                   IceModelVec3 &result) {

  Proj crs(projection, "EPSG:4326");

  IceGrid::ConstPtr grid = result.grid();

  double dx2 = 0.5 * grid->dx(), dy2 = 0.5 * grid->dy();
  double x_offsets[] = {-dx2, dx2, dx2, -dx2};
  double y_offsets[] = {-dy2, -dy2, dy2, dy2};

  IceModelVec::AccessList list{&result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double x0 = grid->x(i), y0 = grid->y(j);

    double *values = result.get_column(i, j);

    for (int k = 0; k < 4; ++k) {

      PJ_COORD in, out;

      in.xy = {x0 + x_offsets[k], y0 + y_offsets[k]};

      // compute lon,lat coordinates:
      out = proj_trans(*crs, PJ_FWD, in);

      if (which == LATITUDE) {
        values[k] = out.lp.lam;
      } else {
        values[k] = out.lp.phi;
      }
    }
  }
}

#else

void compute_cell_areas(const std::string &projection, IceModelVec2S &result) {
  (void) projection;

  IceGrid::ConstPtr grid = result.grid();
  result.set(grid->dx() * grid->dy());
}

static void compute_lon_lat(const std::string &projection, LonLat which,
                            IceModelVec2S &result) {
  (void) projection;
  (void) which;
  (void) result;

  throw RuntimeError(PISM_ERROR_LOCATION, "Cannot compute longitude and latitude."
                     " Please rebuild PISM with PROJ.");
}

static void compute_lon_lat_bounds(const std::string &projection,
                                   LonLat which,
                                   IceModelVec3 &result) {
  (void) projection;
  (void) which;
  (void) result;

  throw RuntimeError(PISM_ERROR_LOCATION, "Cannot compute longitude and latitude bounds."
                     " Please rebuild PISM with PROJ.");
}

#endif

void compute_longitude(const std::string &projection, IceModelVec2S &result) {
  compute_lon_lat(projection, LONGITUDE, result);
}
void compute_latitude(const std::string &projection, IceModelVec2S &result) {
  compute_lon_lat(projection, LATITUDE, result);
}

void compute_lon_bounds(const std::string &projection, IceModelVec3 &result) {
  compute_lon_lat_bounds(projection, LONGITUDE, result);
}

void compute_lat_bounds(const std::string &projection, IceModelVec3 &result) {
  compute_lon_lat_bounds(projection, LATITUDE, result);
}

} // end of namespace pism
