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

#include <cassert>


#include "Lagrange_IO.hh"
#include "base/util/io/PIO.hh"


using namespace pism;


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
void lagrange_prepare_file(const pism::PIO &nc, 
			   unsigned int n_tracers) {
  // FIXME: these should be configurable
  bool dim_exists = false;
  dim_exists = nc.inq_dim("time");
  if (not dim_exists) 
    nc.def_dim("time", pism::PISM_UNLIMITED);

  if (not nc.inq_var("tracer_counter")){
    std::vector<std::string> dims(1);
    dims[0] = "time";
    nc.def_var("tracer_counter", PISM_DOUBLE, dims);
  }

  std::vector<double> index;
  unsigned int count;
  // distributions
  dim_exists = false;
  dim_exists = nc.inq_dim("tracer_index");
  if (dim_exists == false) {
    nc.redef();
    nc.def_dim("tracer_index", n_tracers );
    std::vector<std::string> dims(1);
    dims[0] = "tracer_index";
    nc.def_var("tracer_index", PISM_DOUBLE, dims);
    
    index.resize(n_tracers);
    for (unsigned int j = 0; j < n_tracers; ++j) {
      index[j] = j;
    }
    nc.put_1d_var("tracer_index", 0 , n_tracers, index);
  }


  std::vector<std::string> dims(1);
  dims[0] = "tracer_index";
  // particle positions
  std::vector<std::string> fields ;
  fields.push_back("p_x");
  fields.push_back("p_y");
  fields.push_back("p_z");
  fields.push_back("p_id");
  
  for (std::vector<std::string>::iterator it = fields.begin(); it != fields.end(); ++it) {
    if (not nc.inq_var(*it)) {
      nc.redef();
      nc.def_var(*it, PISM_DOUBLE, dims);
    }
  }
}

void lagrange_prepare_particle_log_file(const pism::PIO &nc){
 
  std::vector<std::string> dims(1);
  dims[0] = "tracer_index";

  bool dim_exists = false;
  dim_exists = nc.inq_dim("tracer_index");
  if (not dim_exists){ 
    nc.def_dim(dims[0], pism::PISM_UNLIMITED); 
    nc.def_var(dims[0], PISM_DOUBLE, dims);
  }
  // particle positions and time -- this is for creation and similar logs

  std::vector<std::string> fields ;
  fields.push_back("p_x");
  fields.push_back("p_y");
  fields.push_back("p_z");
  fields.push_back("p_id");

  
  for (std::vector<std::string>::iterator it = fields.begin(); it != fields.end(); ++it) {
    if (not nc.inq_var(*it)) {
      nc.redef();
      nc.def_var(*it, PISM_DOUBLE, dims);
    }
  }
  
}

