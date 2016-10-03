// Copyright (C) 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"
#include "pismmerge.hh"

#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/error_handling.hh"
#include "base/util/Logger.hh"

using namespace pism;

static char help[] =
  "Tool for merging PISM output files produced using '-o_format quilt'.\n";

using pism::io::NC4_Serial;

int process_one_variable(std::string var_name, std::string input_file, std::string output_file,
                         unsigned int compression_level) {
  NC4_Serial input(MPI_COMM_SELF, 0),
    output(MPI_COMM_SELF, compression_level);
  bool exists;

  fprintf(stderr, "Merging variable %s from %s into %s, compression level %d...\n",
          var_name.c_str(), input_file.c_str(), output_file.c_str(), compression_level);

  // Fill the output file with metadata using the rank=0 "patch".
  input.open(patch_filename(input_file, 0), PISM_READONLY);

  // Create the output file
  output.create(output_file);

  // global attributes
  copy_attributes(input, output, "PISM_GLOBAL");

  define_variable(input, output, var_name);

  input.inq_varid("time_bounds", exists);
  if (exists) {
    define_variable(input, output, "time_bounds");
  }

  // mapping
  input.inq_varid("mapping", exists);
  if (exists) {
    define_variable(input, output, "mapping");
  }

  // pism_override
  input.inq_varid("pism_override", exists);
  if (exists) {
    define_variable(input, output, "pism_override");
  }

  // run_stats
  input.inq_varid("run_stats", exists);
  if (exists) {
    define_variable(input, output, "run_stats");
  }

  // timestamp
  input.inq_varid("timestamp", exists);
  if (exists) {
    define_variable(input, output, "timestamp");
  }

  // lat
  input.inq_varid("lat", exists);
  if (exists) {
    define_variable(input, output, "lat");
  }

  // lon
  input.inq_varid("lon", exists);
  if (exists) {
    define_variable(input, output, "lon");
  }

  input.close();

  copy_all_variables(input_file, output);

  output.close();

  fprintf(stderr, "Done.\n");

  return 0;
}

int process_all_variables(std::string input_file, std::string output_file,
                          unsigned int compression_level) {
  NC4_Serial input(MPI_COMM_SELF, 0),
    output(MPI_COMM_SELF, compression_level);

  fprintf(stderr, "Merging all variables from %s into %s, compression level %d...\n",
          input_file.c_str(), output_file.c_str(), compression_level);

  // Fill the output file with metadata using the rank=0 "patch".
  input.open(patch_filename(input_file, 0), PISM_READONLY);

  // Create the output file
  output.create(output_file);

  // global attributes
  copy_attributes(input, output, "PISM_GLOBAL");

  // define all variables (except for {x,y}_patch)
  int n_vars;
  input.inq_nvars(n_vars);

  for (int j = 0; j < n_vars; ++j) {
    std::string var_name;

    input.inq_varname(j, var_name);

    if (var_name == "x_patch" || var_name == "y_patch") {
      continue;
    }

    define_variable(input, output, var_name);
  }

  copy_all_variables(input_file, output);

  output.close();

  fprintf(stderr, "Done.\n");

  return 0;
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com = MPI_COMM_WORLD;
  int rank;

  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    Logger log(com, 2);

    options::String input_file("-i", "Input file name");
    options::String output_name("-o", "Output file name");
    options::String var_name("-v", "Name of the variable to merge");
    options::Integer compression_level("-L", "Output compression level", 0);
    std::string usage =
      "  Merges output file created using '-o_format quilt'.\n\n"
      "  pismmerge {-i in.nc} [-o out.nc]\n"
      "where:\n"
      "  -i          in.nc is the name used with -extra_file or -o, e.g. ex-RANK.nc\n"
      "  -o          out.nc is the name of the output file that will contain merged data\n"
      "  -v var_name name of the variable to merge\n"
      "  -L <number> output compression level (from 0 to 9)\n"
      "notes:\n"
      "  * -o is optional\n";

    std::vector<std::string> required(1, "-i");

    bool done = show_usage_check_req_opts(log, "PISM-MERGE %s (output file merging tool)",
                                          required, usage);
    if (done) {
      return 0;
    }

    check_input_files(input_file);

    // Check the validity of the -L option.
    if (compression_level.is_set()) {
      if (compression_level < 0 || compression_level > 9) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid compression level: %d.",
                                      compression_level.value());
      }
    }

    // Set the output file name.
    std::string o_name = output_name;
    if (not output_name.is_set()) {
      if (not var_name.is_set()) {
        o_name = output_filename(input_file, "ALL");
      } else {
        o_name = output_filename(input_file, var_name);
      }
    }

    if (rank == 0) {
      if (var_name.is_set()) {
        process_one_variable(var_name, input_file, o_name, compression_level);
      } else {
        process_all_variables(input_file, o_name, compression_level);
      }
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }


  return 0;
}
