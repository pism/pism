/* Copyright (C) 2025 PISM Authors
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

#include <set>

#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/io/OutputWriter.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace io {

void define_variables(const OutputFile &file,
                      const std::set<VariableMetadata> &variables,
                      const VariableMetadata &mapping,
                      bool use_internal_units) {

  std::string mapping_variable_name{};
  if (mapping.has_attributes()) {
    mapping_variable_name = mapping.get_name();
  }

  std::set<std::string> variable_names;
  for (const auto &v : variables) {
    variable_names.insert(v.get_name());
  }

  for (auto var : variables) {
    if (var.get_name() == "PISM_GLOBAL") {
      file.append_history(var["history"]);

      // clear the history attribute
      auto tmp = var;
      tmp["history"] = "";
      file.set_global_attributes(tmp.all_strings(), tmp.all_doubles());

      continue;
    }

    if (use_internal_units) {
      var.output_units(var["units"]);
    }

    auto var_name = var.get_name();

    if (var_name == "lat" and set_member("lat_bnds", variable_names)) {
      var["bounds"] = "lat_bnds";
    }
    if (var_name == "lon" and set_member("lon_bnds", variable_names)) {
      var["bounds"] = "lon_bnds";
    }

    // add extra attributes such as "grid_mapping" and "coordinates". Variables lat, lon,
    // lat_bnds, and lon_bnds should not have these attributes to support CDO (see issue
    // #384).
    //
    // We check names of x and y dimensions to avoid setting extra attributes for variables
    // that use a different grid (e.g. viscous_bed_displacement written by the Lingle-Clark
    // bed deformation model).
    bool have_lat_lon = set_member("lat", variable_names) and set_member("lon", variable_names);
    auto dim_names    = var.dimension_names();
    if (vector_member("x", dim_names) and vector_member("y", dim_names)) {
      if (have_lat_lon and not set_member(var_name, { "lat_bnds", "lon_bnds", "lat", "lon" })) {
        var["coordinates"] = "lat lon";
      }

      if (not mapping_variable_name.empty()) {
        var["grid_mapping"] = mapping_variable_name;
      }
    }

    file.define_variable(var);
  } // end of the loop over variables
}

//! \brief Moves the file aside (file.nc -> file.nc~).
/*!
 * Note: only one processor does the renaming.
 */
void move_if_exists(MPI_Comm com, const std::string &file_to_move, int rank_to_use) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(com, &rank);
  std::string backup_filename = file_to_move + "~";

  if (rank == rank_to_use) {
    bool exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(file_to_move.c_str(), "r")) {
      fclose(f);
      exists = true;
    } else {
      exists = false;
    }

    if (exists) {
      stat = rename(file_to_move.c_str(), backup_filename.c_str());
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, com);

  if (global_stat != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "PISM ERROR: can't move '%s' to '%s'",
                                  file_to_move.c_str(), backup_filename.c_str());
  }
}

//! \brief Check if a file is present are remove it.
/*!
 * Note: only processor 0 does the job.
 */
void remove_if_exists(MPI_Comm com, const std::string &file_to_remove, int rank_to_use) {
  int stat = 0, rank = 0;
  MPI_Comm_rank(com, &rank);

  if (rank == rank_to_use) {
    bool exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(file_to_remove.c_str(), "r")) {
      fclose(f);
      exists = true;
    } else {
      exists = false;
    }

    if (exists) {
      stat = remove(file_to_remove.c_str());
    }
  } // end of "if (rank == rank_to_use)"

  int global_stat = 0;
  MPI_Allreduce(&stat, &global_stat, 1, MPI_INT, MPI_SUM, com);

  if (global_stat != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "PISM ERROR: can't remove '%s'",
                                  file_to_remove.c_str());
  }
}

} // namespace io
} // namespace pism
