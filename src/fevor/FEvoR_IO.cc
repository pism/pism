/* Copyright (C) 2014 PISM Authors
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

#include "FEvoR_IO.hh"
#include "PIO.hh"
#include "NCVariable.hh"

#include "fevor_distribution.hh"

using namespace pism;

int fevor_load_distribution(const std::string &filename,
                            MPI_Comm comm, const UnitSystem &sys,
                            unsigned int index, FEvoR::Distribution &distribution) {
  PetscErrorCode ierr;

  PIO nc(comm, "netcdf3", sys);
  ierr = nc.open(filename, PISM_READWRITE); CHKERRQ(ierr);

  std::vector<double> data(distribution.getNumberCrystals() * FEvoR::Distribution::numberParameters);

  std::string variable_name = "distributions";
  std::vector<unsigned int> start(3), count(3);
  start[0] = index;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = distribution.getNumberCrystals();
  count[2] = FEvoR::Distribution::numberParameters;
  ierr = nc.get_vara_double(variable_name, start, count, &data[0]); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  distribution.loadDistribution(data);

  return 0;
}

/** 
 * Save a distribution to a file at position index.
 *
 * @param filename output file name
 * @param index distribution index in the file
 * @param distribution FEvoR distribution
 *
 * @return 0 on success
 */
int fevor_save_distribution(const std::string &filename,
                            MPI_Comm comm, const UnitSystem &sys,
                            unsigned int index, const FEvoR::Distribution &distribution) {
  PetscErrorCode ierr;

  PIO nc(comm, "netcdf3", sys);
  ierr = nc.open(filename, PISM_READWRITE); CHKERRQ(ierr);

  std::string variable_name = "distributions";
  bool variable_exists = false;
  ierr = nc.inq_var(variable_name, variable_exists); CHKERRQ(ierr);
  if (not variable_exists) {
    std::vector<std::string> dimensions;
    dimensions.push_back("distribution_index");
    dimensions.push_back("crystal_index");
    dimensions.push_back("parameter_index");
    ierr = nc.redef(); CHKERRQ(ierr);
    ierr = nc.def_var(variable_name, PISM_DOUBLE, dimensions); CHKERRQ(ierr);
  }

  std::vector<double> data;
  distribution.saveDistribution(data);

  std::vector<unsigned int> start(3), count(3);
  start[0] = index;
  start[1] = 0;
  start[2] = 0;

  count[0] = 1;
  count[1] = distribution.getNumberCrystals();
  count[2] = FEvoR::Distribution::numberParameters;

  ierr = nc.put_vara_double(variable_name, start, count, &data[0]); CHKERRQ(ierr);

  ierr = nc.close();

  return 0;
}

/**
 * Prepare a file for writing distributions.
 *
 * @param filename output file name
 * @param n_distributions number of distributions
 * @param n_crystals number of crystals per distribution
 *
 * @return 0 on success
 */
int fevor_prepare_file(const std::string &filename,
                       MPI_Comm comm, const UnitSystem &sys,
                       unsigned int n_distributions,
                       unsigned int n_crystals) {
  PetscErrorCode ierr;

  PIO nc(comm, "netcdf3", sys);

  ierr = nc.open(filename, PISM_READWRITE_MOVE); CHKERRQ(ierr);

  std::vector<double> index;
  std::vector<unsigned int> start(1, 0), count(1, 0);
  // distributions
  {
    NCVariable distribution_index("distribution_index", sys, 1);
    ierr = nc.def_dim(n_distributions, distribution_index); CHKERRQ(ierr);

    index.resize(n_distributions);
    for (unsigned int j = 0; j < n_distributions; ++j) {
      index[j] = j;
    }
    count[0] = n_distributions;
    ierr = nc.put_vara_double(distribution_index.get_name(),
                              start, count, &index[0]); CHKERRQ(ierr);
  }

  // crystals
  {
    NCVariable crystal_index("crystal_index", sys, 1);
    ierr = nc.def_dim(n_crystals, crystal_index); CHKERRQ(ierr);

    index.resize(n_crystals);
    for (unsigned int j = 0; j < n_crystals; ++j) {
      index[j] = j;
    }
    count[0] = n_crystals;
    ierr = nc.put_vara_double(crystal_index.get_name(),
                              start, count, &index[0]); CHKERRQ(ierr);
  }

  // parameters
  {
    NCVariable parameter_index("parameter_index", sys, 1);
    ierr = nc.def_dim(FEvoR::Distribution::numberParameters, parameter_index); CHKERRQ(ierr);

    index.resize(FEvoR::Distribution::numberParameters);
    for (unsigned int j = 0; j < FEvoR::Distribution::numberParameters; ++j) {
      index[j] = j;
    }
    count[0] = FEvoR::Distribution::numberParameters;
    ierr = nc.put_vara_double(parameter_index.get_name(),
                              start, count, &index[0]); CHKERRQ(ierr);
  }

  ierr = nc.close();

  return 0;
}
