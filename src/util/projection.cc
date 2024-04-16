/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024 PISM Authors
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
#include <cmath>                // std::pow, std::sqrt, std::fabs
#include <string>
#include <vector>

#include "pism/util/projection.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Array3D.hh"

#include "pism/pism_config.hh"
#include "pism/util/pism_utilities.hh"

#if (Pism_USE_PROJ==1)
#include "pism/util/Proj.hh"
#endif

namespace pism {

MappingInfo::MappingInfo(const std::string &mapping_name, units::System::Ptr unit_system)
  : mapping(mapping_name, unit_system) {
  // empty
}

int parse_epsg(const std::string &proj_string) {
  try {
    int auth_len                    = 5; // length of "epsg:"
    std::string::size_type position = std::string::npos;
    for (const auto &auth : { "epsg:", "EPSG:" }) {
      position = proj_string.find(auth);
      if (position != std::string::npos) {
        break;
      }
    }

    if (position == std::string::npos) {
      throw RuntimeError(PISM_ERROR_LOCATION, "could not find 'epsg:' or 'EPSG:'");
    }


    {
      auto epsg_string = proj_string.substr(position + auth_len);
      char *endptr     = NULL;
      long int result  = strtol(epsg_string.c_str(), &endptr, 10);
      if (endptr == epsg_string.c_str()) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't parse %s (expected an integer)",
                                      epsg_string.c_str());
      }
      return (int)result;
    }

  } catch (RuntimeError &e) {
    e.add_context("parsing the EPSG code in the PROJ string '%s'",
                  proj_string.c_str());
    throw;
  }
}

//! @brief Return CF-Convention "mapping" variable corresponding to an EPSG code specified in a
//! PROJ string.
VariableMetadata epsg_to_cf(units::System::Ptr system, const std::string &proj_string) {
  VariableMetadata mapping("mapping", system);

  int epsg = parse_epsg(proj_string);

  switch (epsg) {
  case 3413:
    mapping["latitude_of_projection_origin"]         = {90.0};
    mapping["scale_factor_at_projection_origin"]     = {1.0};
    mapping["straight_vertical_longitude_from_pole"] = {-45.0};
    mapping["standard_parallel"]                     = {70.0};
    mapping["false_northing"]                        = {0.0};
    mapping["grid_mapping_name"]                     = "polar_stereographic";
    mapping["false_easting"]                         = {0.0};
    break;
  case 3031:
    mapping["latitude_of_projection_origin"]         = {-90.0};
    mapping["scale_factor_at_projection_origin"]     = {1.0};
    mapping["straight_vertical_longitude_from_pole"] = {0.0};
    mapping["standard_parallel"]                     = {-71.0};
    mapping["false_northing"]                        = {0.0};
    mapping["grid_mapping_name"]                     = "polar_stereographic";
    mapping["false_easting"]                         = {0.0};
    break;
  case 3057:
    mapping["grid_mapping_name"]                     = "lambert_conformal_conic" ;
    mapping["longitude_of_central_meridian"]         = {-19.} ;
    mapping["false_easting"]                         = {500000.} ;
    mapping["false_northing"]                        = {500000.} ;
    mapping["latitude_of_projection_origin"]         = {65.} ;
    mapping["standard_parallel"]                     = {64.25, 65.75} ;
    mapping["long_name"]                             = "CRS definition" ;
    mapping["longitude_of_prime_meridian"]           = {0.} ;
    mapping["semi_major_axis"]                       = {6378137.} ;
    mapping["inverse_flattening"]                    = {298.257222101} ;
    break;
  case 5936:
    mapping["latitude_of_projection_origin"]         = {90.0};
    mapping["scale_factor_at_projection_origin"]     = {1.0};
    mapping["straight_vertical_longitude_from_pole"] = {-150.0};
    mapping["standard_parallel"]                     = {90.0};
    mapping["false_northing"]                        = {2000000.0};
    mapping["grid_mapping_name"]                     = "polar_stereographic";
    mapping["false_easting"]                         = {2000000.0};
    break;
  case 26710:
    mapping["longitude_of_central_meridian"]         = {-123.0};
    mapping["false_easting"]                         = {500000.0};
    mapping["false_northing"]                        = {0.0};
    mapping["grid_mapping_name"]                     = "transverse_mercator";
    mapping["inverse_flattening"]                    = {294.978698213898};
    mapping["latitude_of_projection_origin"]         = {0.0};
    mapping["scale_factor_at_central_meridian"]      = {0.9996};
    mapping["semi_major_axis"]                       = {6378206.4};
    mapping["unit"]                                  = "metre";
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
    for (const auto &s : epsg_mapping.all_strings()) {
      if (not info.mapping.has_attribute(s.first)) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistent metadata:\n"
                                      "PROJ string \"%s\" requires %s = \"%s\",\n"
                                      "but the mapping variable has no %s.",
                                      info.proj.c_str(),
                                      s.first.c_str(), s.second.c_str(),
                                      s.first.c_str());
      }

      std::string string = info.mapping[s.first];

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
    for (auto d : epsg_mapping.all_doubles()) {
      if (not info.mapping.has_attribute(d.first)) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has no %s.",
                                      info.proj.c_str(),
                                      d.first.c_str(), d.second[0],
                                      d.first.c_str());
      }

      double value = info.mapping.get_number(d.first);

      if (std::fabs(value - d.second[0]) > 1e-12) {
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

  if (result.proj.empty()) {
    // try the "proj4" attribute
    result.proj = input_file.read_text_attribute("PISM_GLOBAL", "proj4");
  }

  bool proj_is_epsg = false;
  for (const auto &auth : {"epsg:", "EPSG:"}) {
    if (result.proj.find(auth) != std::string::npos) {
      proj_is_epsg = true;
      break;
    }
  }

  if (proj_is_epsg) {
    // re-create the PROJ string to make sure it does not contain "+init="
    result.proj = pism::printf("EPSG:%d", parse_epsg(result.proj));
  }

  if (input_file.variable_exists(mapping_name)) {
    // input file has a mapping variable

    result.mapping = io::read_attributes(input_file, mapping_name, unit_system);

    if (proj_is_epsg) {
      // check consistency
      try {
        check_consistency_epsg(result);
      } catch (RuntimeError &e) {
        e.add_context("getting projection info from %s", input_file.name().c_str());
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
static double triangle_area(const double *A, const double *B, const double *C) {
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

void compute_cell_areas(const std::string &projection, array::Scalar &result) {
  auto grid = result.grid();

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

  array::AccessScope list(result);

  for (auto p = grid->points(); p; p.next()) {
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
                            LonLat which, array::Scalar &result) {

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

  auto grid = result.grid();

  array::AccessScope list{&result};

  for (auto p = grid->points(); p; p.next()) {
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
                                   array::Array3D &result) {

  Proj crs(projection, "EPSG:4326");

  auto grid = result.grid();

  double dx2 = 0.5 * grid->dx(), dy2 = 0.5 * grid->dy();
  double x_offsets[] = {-dx2, dx2, dx2, -dx2};
  double y_offsets[] = {-dy2, -dy2, dy2, dy2};

  array::AccessScope list{&result};

  for (auto p = grid->points(); p; p.next()) {
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

void compute_cell_areas(const std::string &projection, array::Scalar &result) {
  (void) projection;

  auto grid = result.grid();
  result.set(grid->dx() * grid->dy());
}

static void compute_lon_lat(const std::string &projection, LonLat which,
                            array::Scalar &result) {
  (void) projection;
  (void) which;
  (void) result;

  throw RuntimeError(PISM_ERROR_LOCATION, "Cannot compute longitude and latitude."
                     " Please rebuild PISM with PROJ.");
}

static void compute_lon_lat_bounds(const std::string &projection,
                                   LonLat which,
                                   array::Array3D &result) {
  (void) projection;
  (void) which;
  (void) result;

  throw RuntimeError(PISM_ERROR_LOCATION, "Cannot compute longitude and latitude bounds."
                     " Please rebuild PISM with PROJ.");
}

#endif

void compute_longitude(const std::string &projection, array::Scalar &result) {
  compute_lon_lat(projection, LONGITUDE, result);
}
void compute_latitude(const std::string &projection, array::Scalar &result) {
  compute_lon_lat(projection, LATITUDE, result);
}

void compute_lon_bounds(const std::string &projection, array::Array3D &result) {
  compute_lon_lat_bounds(projection, LONGITUDE, result);
}

void compute_lat_bounds(const std::string &projection, array::Array3D &result) {
  compute_lon_lat_bounds(projection, LATITUDE, result);
}

static std::string albers_conical_equal_area_to_proj(const VariableMetadata &mapping) {
  // standard_parallel - There may be 1 or 2 values.
  // longitude_of_central_meridian
  // latitude_of_projection_origin
  // false_easting
  // false_northing
  double longitude_of_projection_origin = mapping["longitude_of_projection_origin"];
  double latitude_of_projection_origin  = mapping["latitude_of_projection_origin"];
  double false_easting                  = mapping["false_easting"];
  double false_northing                 = mapping["false_northing"];

  std::vector<double> standard_parallel = mapping["standard_parallel"];

  auto result = pism::printf("+proj=aea +type=crs +lon_0=%f +lat_0=%f +x_0=%f +y_0=%f +lat_1=%f",
                             longitude_of_projection_origin,
                             latitude_of_projection_origin,
                             false_easting,
                             false_northing,
                             standard_parallel[0]);
  if (standard_parallel.size() == 2) {
    result += pism::printf(" +lat_2=%f", standard_parallel[1]);
  }

  return result;
}

static std::string azimuthal_equidistant_to_proj(const VariableMetadata &mapping) {
  // Parameters:
  // longitude_of_projection_origin
  // latitude_of_projection_origin
  // false_easting
  // false_northing
  double longitude_of_projection_origin = mapping["longitude_of_projection_origin"];
  double latitude_of_projection_origin  = mapping["latitude_of_projection_origin"];
  double false_easting                  = mapping["false_easting"];
  double false_northing                 = mapping["false_northing"];
  return pism::printf("+proj=aeqd +type=crs +lon_0=%f +lat_0=%f +x_0=%f +y_0=%f",
                      longitude_of_projection_origin,
                      latitude_of_projection_origin, false_easting, false_northing);
}

static std::string lambert_azimuthal_equal_area_to_proj(const VariableMetadata &mapping) {
  // longitude_of_projection_origin
  // latitude_of_projection_origin
  // false_easting
  // false_northing
  double longitude_of_projection_origin = mapping["longitude_of_projection_origin"];
  double latitude_of_projection_origin  = mapping["latitude_of_projection_origin"];
  double false_easting                  = mapping["false_easting"];
  double false_northing                 = mapping["false_northing"];
  return pism::printf("+proj=laea +type=crs +lon_0=%f +lat_0=%f +x_0=%f +y_0=%f",
                      longitude_of_projection_origin,
                      latitude_of_projection_origin, false_easting, false_northing);
}

static std::string lambert_conformal_conic_to_proj(const VariableMetadata &mapping) {
  // standard_parallel - There may be 1 or 2 values.
  // longitude_of_central_meridian
  // latitude_of_projection_origin
  // false_easting
  // false_northing
  double longitude_of_projection_origin = mapping["longitude_of_projection_origin"];
  double latitude_of_projection_origin  = mapping["latitude_of_projection_origin"];
  double false_easting                  = mapping["false_easting"];
  double false_northing                 = mapping["false_northing"];
  std::vector<double> standard_parallel = mapping["standard_parallel"];

  auto result = pism::printf("+proj=lcc +type=crs +lon_0=%f +lat_0=%f +x_0=%f +y_0=%f +lat_1=%f",
                             longitude_of_projection_origin,
                             latitude_of_projection_origin,
                             false_easting,
                             false_northing,
                             standard_parallel[0]);
  if (standard_parallel.size() == 2) {
    result += pism::printf(" +lat_2=%f", standard_parallel[1]);
  }
  return result;
}

static std::string lambert_cylindrical_equal_area_to_proj(const VariableMetadata &mapping) {
  // Parameters:
  // longitude_of_central_meridian
  // Either standard_parallel or scale_factor_at_projection_origin
  // false_easting
  // false_northing

  double longitude_of_central_meridian = mapping["longitude_of_central_meridian"];
  double false_easting                 = mapping["false_easting"];
  double false_northing                = mapping["false_northing"];

  auto result = pism::printf("+proj=cea +type=crs +lon_0=%f +x_0=%f +y_0=%f",
                             longitude_of_central_meridian, false_easting, false_northing);

  if (mapping.has_attribute("standard_parallel")) {
    result += pism::printf(" +lat_ts=%f", (double)mapping["standard_parallel"]);
  } else {
    result += pism::printf(" +k_0=%f", (double)mapping["scale_factor_at_projection_origin"]);
  }

  return result;
}

static std::string latitude_longitude_to_proj(const VariableMetadata &mapping) {
  // No parameters
  return "+proj=lonlat +type=crs";
}

static std::string mercator_to_proj(const VariableMetadata &mapping) {
  // Parameters:
  // longitude_of_projection_origin
  // Either standard_parallel (EPSG 9805) or scale_factor_at_projection_origin (EPSG 9804)
  // false_easting
  // false_northing

  double longitude_of_projection_origin = mapping["longitude_of_projection_origin"];
  double false_easting                  = mapping["false_easting"];
  double false_northing                 = mapping["false_northing"];

  auto result = pism::printf("+proj=merc +type=crs +lon_0=%f +x_0=%f +y_0=%f",
                             longitude_of_projection_origin,
                             false_easting, false_northing);

  if (mapping.has_attribute("standard_parallel")) {
    result += pism::printf(" +lat_ts=%f", (double)mapping["standard_parallel"]);
  } else {
    result += pism::printf(" +k_0=%f", (double)mapping["scale_factor_at_projection_origin"]);
  }

  return result;
}

static std::string orthographic_to_proj(const VariableMetadata &mapping) {
  // Parameters:
  // longitude_of_projection_origin
  // latitude_of_projection_origin
  // false_easting
  // false_northing
  double longitude_of_projection_origin = mapping["longitude_of_projection_origin"];
  double latitude_of_projection_origin  = mapping["latitude_of_projection_origin"];
  double false_easting                  = mapping["false_easting"];
  double false_northing                 = mapping["false_northing"];

  // clang-format off
  return pism::printf("+proj=ortho +type=crs +lon_0=%f +lat_0=%f +x_0=%f +y_0=%f",
                      longitude_of_projection_origin,
                      latitude_of_projection_origin,
                      false_easting,
                      false_northing);
  // clang-format on
}

static std::string polar_stereographic_to_proj(const VariableMetadata &mapping) {
  /* Parameters:
    straight_vertical_longitude_from_pole
    latitude_of_projection_origin - Either +90. or -90.
    Either standard_parallel or scale_factor_at_projection_origin
    false_easting
    false_northing
  */

  double straight_vertical_longitude_from_pole = mapping["straight_vertical_longitude_from_pole"];
  double latitude_of_projection_origin         = mapping["latitude_of_projection_origin"];
  double false_easting                         = mapping["false_easting"];
  double false_northing                        = mapping["false_northing"];

  auto result = pism::printf("+proj=stere +type=crs +lat_0=%f +lon_0=%f +x_0=%f +y_0=%f",
                             latitude_of_projection_origin,
                             straight_vertical_longitude_from_pole,
                             false_easting,
                             false_northing);

  if (mapping.has_attribute("standard_parallel")) {
    result += pism::printf(" +lat_ts=%f", (double)mapping["standard_parallel"]);
  } else {
    result += pism::printf(" +k_0=%f", (double)mapping["scale_factor_at_projection_origin"]);
  }

  return result;
}

static std::string rotated_latitude_longitude_to_proj(const VariableMetadata &mapping) {
  // grid_north_pole_latitude
  // grid_north_pole_longitude
  // north_pole_grid_longitude - This parameter is optional (default is 0).

  /*
   * From https://github.com/OSGeo/PROJ/blob/356e255022cea47bf242bc9e6345a1e87ab1031d/src/iso19111/operation/conversion.cpp#L2419
   *
   * Several conventions for the pole rotation method exists.
   * The parameters provided in this method are remapped to the PROJ ob_tran
   * operation with:
   * <pre>
   * +proj=ob_tran +o_proj=longlat +o_lon_p=northPoleGridLongitude
   *                               +o_lat_p=gridNorthPoleLatitude
   *                               +lon_0=180+gridNorthPoleLongitude
   * </pre>
   */

  double grid_north_pole_latitude  = mapping["grid_north_pole_latitude"];
  double grid_north_pole_longitude = mapping["grid_north_pole_longitude"];
  double north_pole_grid_longitude = 0.0;

  if (mapping.has_attribute("north_pole_grid_longitude")) {
    north_pole_grid_longitude = mapping["north_pole_grid_longitude"];
  }

  return pism::printf("+proj=ob_tran +type=crs +o_proj=longlat +o_lon_p=%f +o_lat_p=%f +lon_0=%f",
                      north_pole_grid_longitude,
                      grid_north_pole_latitude,
                      180.0 + grid_north_pole_longitude);
}

static std::string stereographic_to_proj(const VariableMetadata &mapping) {
  /* Parameters:
    longitude_of_projection_origin
    latitude_of_projection_origin
    scale_factor_at_projection_origin
    false_easting
    false_northing
*/

  double straight_vertical_longitude_from_pole = mapping["straight_vertical_longitude_from_pole"];
  double latitude_of_projection_origin         = mapping["latitude_of_projection_origin"];
  double false_easting                         = mapping["false_easting"];
  double false_northing                        = mapping["false_northing"];
  double scale_factor_at_projection_origin     = mapping["scale_factor_at_projection_origin"];

  // clang-format off
  return pism::printf("+proj=stere +type=crs +lat_0=%f +lon_0=%f +x_0=%f +y_0=%f +k_0=%f",
                      latitude_of_projection_origin,
                      straight_vertical_longitude_from_pole,
                      false_easting,
                      false_northing,
                      scale_factor_at_projection_origin);
  // clang-format on
}

static std::string transverse_mercator_to_proj(const VariableMetadata &mapping) {
  double scale_factor_at_central_meridian = mapping["scale_factor_at_central_meridian"];
  double longitude_of_central_meridian    = mapping["longitude_of_central_meridian"];
  double latitude_of_projection_origin    = mapping["latitude_of_projection_origin"];
  double false_easting                    = mapping["false_easting"];
  double false_northing                   = mapping["false_northing"];

  return pism::printf("+proj=tmerc +type=crs +lon_0=%f +lat_0=%f +k=%f +x_0=%f +y_0=%f",
                      longitude_of_central_meridian,
                      latitude_of_projection_origin,
                      scale_factor_at_central_meridian,
                      false_easting,
                      false_northing);
}

static std::string vertical_perspective_to_proj(const VariableMetadata &mapping) {
  // Parameters:
  // latitude_of_projection_origin
  // longitude_of_projection_origin
  // perspective_point_height
  // false_easting
  // false_northing
  std::string name = mapping["grid_mapping_name"];
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "grid mapping '%s' is not supported",
                                name.c_str());
  return "";
}

static std::string common_flags(const VariableMetadata &mapping) {

  auto flag = [&mapping](const char* name, const char* proj_name) {
    if (mapping.has_attribute(name)) {
      return pism::printf(" +%s=%f", proj_name, (double)mapping[name]);
    }
    return std::string{};
  };

  std::string result;

  // ellipsoid
  result += flag("earth_radius", "R");
  result += flag("inverse_flattening", "rf");
  result += flag("semi_major_axis", "a");
  result += flag("semi_minor_axis", "b");

  // prime meridian
  result += flag("longitude_of_prime_meridian", "pm");

  return result;
}

/*!
 * Convert a CF-style grid mapping definition to a PROJ string.
 *
 * See http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/apf.html
 *
 * See also https://proj.org/en/9.4/operations/projections/index.html for a list of
 * projections and their parameters.
 *
 * PROJ strings returned by this function are inspired by
 * https://scitools.org.uk/cartopy/docs/latest/index.html
 */
std::string cf_to_proj(const VariableMetadata &mapping) {

  // all projections: use the crs_wkt attribute if it is present
  std::string wkt = mapping["crs_wkt"];
  if (not wkt.empty()) {
    return wkt;
  }

  typedef std::string(*cf_to_proj_fn)(const VariableMetadata&);

  std::map<std::string, cf_to_proj_fn> mappings = {
    { "albers_conical_equal_area", albers_conical_equal_area_to_proj },
    { "azimuthal_equidistant", azimuthal_equidistant_to_proj },
    { "lambert_azimuthal_equal_area", lambert_azimuthal_equal_area_to_proj },
    { "lambert_conformal_conic", lambert_conformal_conic_to_proj },
    { "lambert_cylindrical_equal_area", lambert_cylindrical_equal_area_to_proj },
    { "latitude_longitude", latitude_longitude_to_proj },
    { "mercator", mercator_to_proj },
    { "orthographic", orthographic_to_proj },
    { "polar_stereographic", polar_stereographic_to_proj },
    { "rotated_latitude_longitude", rotated_latitude_longitude_to_proj },
    { "stereographic", stereographic_to_proj },
    { "transverse_mercator", transverse_mercator_to_proj },
    { "vertical_perspective", vertical_perspective_to_proj },
  };

  std::string mapping_name = mapping["grid_mapping_name"];

  if (mappings.find(mapping_name) == mappings.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "conversion CF -> PROJ for a '%s' grid is not implemented",
                                  mapping_name.c_str());
  }

  return mappings[mapping_name](mapping) + common_flags(mapping);
}



} // end of namespace pism
