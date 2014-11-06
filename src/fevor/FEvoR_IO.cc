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

int fevor_load_distribution(const PIO &nc,
                            unsigned int index,
                            unsigned int time_index,
                            FEvoR::Distribution &distribution) {
  PetscErrorCode ierr;

  size_t data_size = distribution.getNumberCrystals() * FEvoR::Distribution::numberParameters;
  std::vector<double> data(data_size);

  std::vector<unsigned int> start(4), count(4);
  start[0] = time_index;
  start[1] = index;
  start[2] = 0;
  start[3] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = distribution.getNumberCrystals();
  count[3] = FEvoR::Distribution::numberParameters;

  ierr = nc.get_vara_double("distributions", start, count, &data[0]); CHKERRQ(ierr);

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
int fevor_save_distribution(const PIO &nc,
                            unsigned int index,
                            unsigned int time_index,
                            const FEvoR::Distribution &distribution) {
  PetscErrorCode ierr;

  std::vector<double> data;
  distribution.saveDistribution(data);

  std::vector<unsigned int> start(4), count(4);
  start[0] = time_index;
  start[1] = index;
  start[2] = 0;
  start[3] = 0;

  count[0] = 1;
  count[1] = 1;
  count[2] = distribution.getNumberCrystals();
  count[3] = FEvoR::Distribution::numberParameters;

  ierr = nc.put_vara_double("distributions", start, count, &data[0]); CHKERRQ(ierr);

  return 0;
}

/**
 * Prepare a file for writing distributions.
 *
 * @param filename output file name
 * @param n_distributions number of distributions
 * @param packing_dimensions number of crystals per distribution in each spatial direction
 *
 * Total number of crystals per distribution is equal to the product
 * of numbers in packing_dimensions.
 *
 * @return 0 on success
 */
int fevor_prepare_file(const pism::PIO &nc, const pism::UnitSystem &sys,
                       unsigned int n_distributions,
                       const std::vector<unsigned int> packing_dimensions) {
  PetscErrorCode ierr;

  // FIXME: these should be configurable
  ierr = nc.def_time("time", "standard", "seconds"); CHKERRQ(ierr);

  // save packing dimensions
  {
    // make a copy that uses "double"
    std::vector<double> pd(packing_dimensions.size());
    for (unsigned int i = 0; i < pd.size(); ++i) {
      pd[i] = packing_dimensions[i];
    }
    // save packing dimensions
    ierr = nc.put_att_double("PISM_GLOBAL", "packing_dimensions",
                             PISM_INT, pd); CHKERRQ(ierr);
  }

  assert(packing_dimensions.size() == 3);
  unsigned int n_crystals = packing_dimensions[0] * packing_dimensions[1] * packing_dimensions[2];

  std::vector<double> index;
  std::vector<unsigned int> start(1, 0), count(1, 0);
  // distributions
  bool dim_exists = false;
  ierr = nc.inq_dim("distribution_index", dim_exists); CHKERRQ(ierr);
  if (dim_exists == false) {
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
  ierr = nc.inq_dim("crystal_index", dim_exists); CHKERRQ(ierr);
  if (dim_exists == false) {
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
  ierr = nc.inq_dim("parameter_index", dim_exists); CHKERRQ(ierr);
  if (dim_exists == false)  {
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

  {
    std::string variable_name = "distributions";
    bool variable_exists = false;
    ierr = nc.inq_var(variable_name, variable_exists); CHKERRQ(ierr);
    if (not variable_exists) {
      std::vector<std::string> dimensions;
      dimensions.push_back("time");
      dimensions.push_back("distribution_index");
      dimensions.push_back("crystal_index");
      dimensions.push_back("parameter_index");
      ierr = nc.redef(); CHKERRQ(ierr);
      ierr = nc.def_var(variable_name, PISM_DOUBLE, dimensions); CHKERRQ(ierr);
    }
  }

  std::vector<std::string> dims(2);
  dims[0] = "time";
  dims[1] = "distribution_index";
  // particle positions
  {
    bool exists = false;
    ierr = nc.inq_var("p_x", exists); CHKERRQ(ierr);
    if (not exists) {
      ierr = nc.def_var("p_x", PISM_DOUBLE, dims); CHKERRQ(ierr);
    }

    ierr = nc.inq_var("p_y", exists); CHKERRQ(ierr);
    if (not exists) {
      ierr = nc.def_var("p_y", PISM_DOUBLE, dims); CHKERRQ(ierr);
    }

    ierr = nc.inq_var("p_z", exists); CHKERRQ(ierr);
    if (not exists) {
      ierr = nc.def_var("p_z", PISM_DOUBLE, dims); CHKERRQ(ierr);
    }
  }

  // diagnostics
  {
    bool exists = false;
    ierr = nc.inq_var("n_migration_recrystallizations", exists); CHKERRQ(ierr);
    if (not exists) {
      ierr = nc.def_var("n_migration_recrystallizations", PISM_DOUBLE, dims); CHKERRQ(ierr);
    }
    ierr = nc.inq_var("n_polygonization_recrystallizations", exists); CHKERRQ(ierr);
    if (not exists) {
      ierr = nc.def_var("n_polygonization_recrystallizations", PISM_DOUBLE, dims); CHKERRQ(ierr);
    }
  }

  return 0;
}

int fevor_load_particle_positions(const pism::PIO &nc,
                                  unsigned int time_index,
                                  std::vector<double> &x,
                                  std::vector<double> &y,
                                  std::vector<double> &z) {
  PetscErrorCode ierr;
  assert(x.size() == y.size());
  assert(x.size() == z.size());

  std::vector<unsigned int> start(2), count(2);
  start[0] = time_index;
  start[1] = 0;

  count[0] = 1;
  count[1] = x.size();

  ierr = nc.get_vara_double("p_x", start, count, &x[0]); CHKERRQ(ierr);
  ierr = nc.get_vara_double("p_y", start, count, &y[0]); CHKERRQ(ierr);
  ierr = nc.get_vara_double("p_z", start, count, &z[0]); CHKERRQ(ierr);

  return 0;
}

int fevor_save_particle_positions(const PIO &nc,
                                  unsigned int time_index,
                                  std::vector<double> &x,
                                  std::vector<double> &y,
                                  std::vector<double> &z) {
  PetscErrorCode ierr;
  assert(x.size() == y.size());
  assert(x.size() == z.size());

  std::vector<unsigned int> start(2), count(2);
  start[0] = time_index;
  start[1] = 0;

  count[0] = 1;
  count[1] = x.size();

  ierr = nc.put_vara_double("p_x", start, count, &x[0]); CHKERRQ(ierr);
  ierr = nc.put_vara_double("p_y", start, count, &y[0]); CHKERRQ(ierr);
  ierr = nc.put_vara_double("p_z", start, count, &z[0]); CHKERRQ(ierr);

  return 0;
}


int fevor_save_recrystallization_numbers(const PIO &nc,
                                         unsigned int time_index,
                                         std::vector<double> &migration,
                                         std::vector<double> &polygonization) {
  PetscErrorCode ierr;
  assert(migration.size() == polygonization.size());

  std::vector<unsigned int> start(2), count(2);
  start[0] = time_index;
  start[1] = 0;

  count[0] = 1;
  count[1] = migration.size();

  ierr = nc.put_vara_double("n_migration_recrystallizations", start, count,
                            &migration[0]); CHKERRQ(ierr);

  ierr = nc.put_vara_double("n_polygonization_recrystallizations", start, count,
                            &polygonization[0]); CHKERRQ(ierr);
  return 0;
}
