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
#include <vector>
#include <algorithm>

#include <fevor_distribution.hh>
#include <vector_tensor_operations.hh>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include "PISMFEvoR.hh"
#include "PISMConfig.hh"
#include "PISMVars.hh"
#include "PISMStressBalance_diagnostics.hh"
#include "enthalpyConverter.hh"

#include "FEvoR_IO.hh"
#include "pism_options.hh"

namespace pism {

PISMFEvoR::PISMFEvoR(IceGrid &g, const Config &conf, EnthalpyConverter *EC,
                     StressBalance *stress_balance)
  : Component_TS(g, conf), m_stress_balance(stress_balance), m_EC(EC),
    m_packing_dimensions(std::vector<unsigned int>(3, 5)),
    m_d_iso(m_packing_dimensions, 0.0)
{

  assert(m_EC != NULL);
  assert(m_stress_balance != NULL);

  PetscErrorCode ierr = allocate(); CHKERRCONTINUE(ierr);

  // will be allocated in init()
  m_pressure = NULL;
  m_tauxz    = NULL;
  m_tauyz    = NULL;

  // points to storage owned by IceModel
  m_enthalpy = NULL;
}

PISMFEvoR::~PISMFEvoR() {
  delete m_pressure;
  delete m_tauxz;
  delete m_tauyz;
}

PetscErrorCode PISMFEvoR::max_timestep(double t, double &dt, bool &restrict) {
  // FIXME: put real code here
  PetscErrorCode ierr = Component_TS::max_timestep(t, dt, restrict); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMFEvoR::update(double t, double dt) {
  m_t = t;
  m_dt = dt;
  PetscErrorCode ierr;

  IceModelVec3 *u = NULL, *v = NULL, *w = NULL;
  ierr = m_stress_balance->get_3d_velocity(u, v, w); CHKERRQ(ierr);

  assert(m_pressure != NULL);
  IceModelVec* pressure = NULL; // FIXME: use a smart pointer
  ierr = m_pressure->compute(pressure); CHKERRQ(ierr);
  IceModelVec3 *pressure3 = static_cast<IceModelVec3*>(pressure);

  assert(m_tauxz != NULL);
  IceModelVec* tauxz = NULL;    // FIXME: use a smart pointer
  ierr = m_tauxz->compute(tauxz); CHKERRQ(ierr);
  IceModelVec3 *tauxz3 = static_cast<IceModelVec3*>(tauxz);

  assert(m_tauyz != NULL);
  IceModelVec* tauyz = NULL;    // FIXME: use a smart pointer
  ierr = m_tauyz->compute(tauyz); CHKERRQ(ierr);
  IceModelVec3 *tauyz3 = static_cast<IceModelVec3*>(tauyz);

  {
    IceModelVec::AccessList list;
    list.add(*u);
    list.add(*v);
    list.add(*w);
    list.add(*pressure3);
    list.add(*tauxz3);
    list.add(*tauyz3);
    list.add(*m_enthalpy);

    unsigned int n_particles = m_distributions.size();

    // Ensure that we have enough storage for diagnostics
    m_n_migration_recrystallizations.resize(n_particles);
    m_n_polygonization_recrystallizations.resize(n_particles);

    // Get enhancement factor for every particle
    for (unsigned int i = 0; i < n_particles; ++i) {

      // interpolate these values from PISM
      double P = 0.0,
        txz    = 0.0,
        tyz    = 0.0,
        E      = 0.0,
        T      = 0.0;

      // check if the point (x,y,z) is within the domain:
      assert(0.0 <= m_p_z[i] && m_p_z[i] <= grid.Lz);
      assert(grid.x.front() <= m_p_x[i] && m_p_x[i] <= grid.x.back());
      assert(grid.y.front() <= m_p_y[i] && m_p_y[i] <= grid.y.back());

      ierr = evaluate_at_point(*pressure3, m_p_x[i], m_p_y[i], m_p_z[i], P); CHKERRQ(ierr);
      ierr = evaluate_at_point(*tauxz3, m_p_x[i], m_p_y[i], m_p_z[i], txz);  CHKERRQ(ierr);
      ierr = evaluate_at_point(*tauyz3, m_p_x[i], m_p_y[i], m_p_z[i], tyz);  CHKERRQ(ierr);

      /* Indexing: {0, 1, 2,
       *            3, 4, 5,
       *            6, 7, 8}
       */
      std::vector<double> stress(9,0.0);
      stress[0] = stress[4] = stress[8] = -P;
      stress[2] = stress[6] = txz;           
      stress[5] = stress[7] = tyz;           

      ierr = evaluate_at_point(*m_enthalpy, m_p_x[i], m_p_y[i], m_p_z[i], E); CHKERRQ(ierr);
      ierr = m_EC->getAbsTemp(E, P, T); CHKERRQ(ierr);

      std::vector<double> bulkEdot(9, 0.0);
      // http://en.wikipedia.org/wiki/Step_in_Time
      m_distributions[i].stepInTime(T, stress, m_t, m_dt,
                                    m_n_migration_recrystallizations[i],
                                    m_n_polygonization_recrystallizations[i],
                                    bulkEdot);

      std::vector<double> bulkEdot_iso(9, 0.0);
      m_d_iso.stepInTime(T, stress, m_t, m_dt, bulkEdot_iso);

      const double M = FEvoR::tensorMagnitude(bulkEdot),
        M_iso = FEvoR::tensorMagnitude(bulkEdot_iso);

      if (M_iso != 0.0) {
        m_p_e[i] = M / M_iso;
      } else {
        // If a particle is outside the ice blob, then pressure == 0
        // and M_iso == 0, so we need this special case.
        m_p_e[i] = 1.0;
      }

      // some bounds for the enhancement factor
      if (m_p_e[i] < 0.4) {
        m_p_e[i] = 0.4;
      } else if (m_p_e[i] > 10.0) {
        m_p_e[i] = 10.0;
        // upper bound.
      }
    } // end of the for-loop over particles

    // set the enhancement factor for every grid point from our particle cloud
    ierr = pointcloud_to_grid(m_p_x, m_p_z, m_p_e, m_enhancement_factor); CHKERRQ(ierr);

    ierr = m_enhancement_factor.update_ghosts(); CHKERRQ(ierr);

    // update particle positions -- don't want to update until gridded
    // enhancement factor is updated
    for (unsigned int i = 0; i < n_particles; ++i) {
      double p_u = 0.0,
             p_v = 0.0,
             p_w = 0.0;

      ierr = evaluate_at_point(*u, m_p_x[i], m_p_y[i], m_p_z[i], p_u);   CHKERRQ(ierr);
      ierr = evaluate_at_point(*v, m_p_x[i], m_p_y[i], m_p_z[i], p_v);   CHKERRQ(ierr);
      ierr = evaluate_at_point(*w, m_p_x[i], m_p_y[i], m_p_z[i], p_w);   CHKERRQ(ierr);

      ierr = update_particle_position(m_p_x[i], m_p_y[i], m_p_z[i], p_u, p_v, p_w, m_dt); CHKERRQ(ierr);
    }
  }

  delete pressure;
  delete tauxz;
  delete tauyz;

  return 0;
}

/**
 * Evaluate new particle position from velocity
 *
 * @param x, y, z coordinates of a point within the domain
 * @param u, v, w velocities of x, y, z
 * @param dt time step length
 *
 * @return 0 on success
 */
PetscErrorCode PISMFEvoR::update_particle_position(double &x, double &y, double &z,
                                                   double u, double v, double w,
                                                   double dt) {

  // basic Euler method. Assuming u, v, w are 'good'
  x += u*dt;
  y += v*dt; // probably not needed and v should be zero in 2D flow line model
  z += w*dt; // w should be zero

  if (z > grid.Lz) {
    z = grid.Lz; 
  } else if (z < 0.0) {
    z = 0.0;
  }
  
  // assuming periodic x... 
  if (x < grid.x.front()) {
    x += grid.x.back()-grid.x.front();
  } else if (x > grid.x.back()) {
    x -= grid.x.back()-grid.x.front();
  }
  
  // assuming periodic y... 
  if (y < grid.y.front()) {
    y += grid.y.back()-grid.y.front();
  } else if (y > grid.y.back()) {
    y -= grid.y.back()-grid.y.front();
  }

  // check if the point (x,y,z) is within the domain:
  assert(0.0 <= z && z <= grid.Lz);
  // PISM's horizontal grid is cell-centered, so (-grid.Lx,
  // -grid.Ly) is actually outside of the convex hull of all grid
  // points.
  assert(grid.x.front() <= x && x <= grid.x.back());
  assert(grid.y.front() <= y && y <= grid.y.back());

  return 0;
}

/**
 * Evaluate a 3D field at a given point.
 *
 * @param input a 3D field
 * @param x, y, z coordinates of a point within the domain
 * @param result interpolated value
 *
 * @return 0 on success
 */
PetscErrorCode PISMFEvoR::evaluate_at_point(IceModelVec3 &input,
                                            double x, double y, double z,
                                            double &result) {
  PetscErrorCode ierr;

  // compute indexes for accessing neighboring columns:
  int I_left = 0, I_right = 0, J_bottom = 0, J_top = 0;
  grid.compute_point_neighbors(x, y, I_left, I_right, J_bottom, J_top);

  // compute index of the vertical level just below z:
  unsigned int K = 0;
  // K + 1 (used below) should be at most Mz - 1
  while (K + 1 < grid.Mz - 1 && grid.zlevels[K + 1] < z) {
    K++;
  }

  // get pointers to neighboring columns:
  double* column[4] = {NULL, NULL, NULL, NULL};
  ierr = input.getInternalColumn(I_left,  J_bottom, &column[0]); CHKERRQ(ierr);
  ierr = input.getInternalColumn(I_right, J_bottom, &column[1]); CHKERRQ(ierr);
  ierr = input.getInternalColumn(I_right, J_top,    &column[2]); CHKERRQ(ierr);
  ierr = input.getInternalColumn(I_left,  J_top,    &column[3]); CHKERRQ(ierr);

  // compute interpolation weights
  std::vector<double> weights = grid.compute_interp_weights(x, y);
  assert(K < grid.Mz);
  const double z_weight = (z - grid.zlevels[K]) / (grid.zlevels[K+1] - grid.zlevels[K]);

  // interpolate
  result = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    // vertical interpolation:
    double f_i = column[i][K] + z_weight * (column[i][K+1] - column[i][K]);
    // horizontal interpolation:
    result += weights[i] * f_i;
  }

  return 0;
}

// typedefs for pointcloud_to_grid()
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>   Traits;
// type to represent 2D points
typedef K::Point_2 Point;
// how to compare our 2D points
typedef K::Less_xy_2 Map_compare;
// field number type -- has models for what we construct our points out of
typedef K::FT Field_type;
// finding the neighbors
typedef std::vector< std::pair< Point, Field_type  > > Point_coordinate_vector;
typedef CGAL::Triple< std::back_insert_iterator<Point_coordinate_vector>,Field_type, bool > The_neighbors;

/**
 * Interpolate values defined at locations (x, z) to the grid and store in result.
 *
 * @param x, z arrays of coordinates of points in the cloud.
 * @param values values at locations (x, z)
 * @param result output field
 *
 * @return 0 on success
 */
PetscErrorCode PISMFEvoR::pointcloud_to_grid(const std::vector<double> &x,
                                             const std::vector<double> &z,
                                             const std::vector<double> &values,
                                             IceModelVec3 &result) {

  assert(x.size() ==z.size() && z.size() == values.size());

  const unsigned int n_particles = x.size();

  Delaunay_triangulation D_TRI;

  // map our points to our function values
  std::map<Point, Field_type, Map_compare> function_values;
  // function to access our data
  typedef CGAL::Data_access< std::map<Point, Field_type, Map_compare > > Value_access;

  // get the points for our convex hull
  for (unsigned int pn = 0; pn < n_particles; ++pn) {
    Point p(x[pn],z[pn]);
    D_TRI.insert( p );
    function_values.insert( std::make_pair(p, values[pn]) );
  }

  IceModelVec::AccessList list(result);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (unsigned int k=0; k<grid.Mz; ++k) {
      Point INTERP( grid.x[i],grid.zlevels[k] );

      // make a vector of the iterpolation point and type
      Point_coordinate_vector coord;
      The_neighbors norm = CGAL::natural_neighbor_coordinates_2 (D_TRI, INTERP, std::back_inserter(coord) );
      
      Field_type res;
      // make sure point is in convex hull!
      if(!norm.third) {
        // FIXME this is dumb logic. Works for only the Enhancement factor now!!
        if (grid.zlevels[k] > grid.zlevels[grid.Mz-1]/2.0) {
          // above the middle
          res = Field_type(1.0); // isotropic.
        } else {
          res = Field_type( *std::max_element(values.begin(), values.end()) ); // max.
        } 
      } else {
        res =  CGAL::linear_interpolation (coord.begin(), coord.end(), norm.second, Value_access(function_values));
      }
      
      // set result for all y grid points at INTERP(x,z)
      for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
        result(i,j,k) = double(res);
      }
    }
  }
  
  // deal with side boundaries. 
  for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
    for (unsigned int k=0; k<grid.Mz; ++k) {
      result(grid.xs,j,k) = result(grid.xs+1,j,k);
      result(grid.xs+grid.xm-1,j,k) = result(grid.xs+grid.xm-2,j,k);
    }
  } 
  return 0;
}

void PISMFEvoR::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword != "none") {
    result.insert(m_enhancement_factor.metadata().get_string("short_name"));
    result.insert("distributions");
  }
}

PetscErrorCode PISMFEvoR::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                           IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "enhancement_factor")) {
    ierr = m_enhancement_factor.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "distributions") or set_contains(vars, "recrystallizations")) {
    ierr = fevor_prepare_file(nc, grid.get_unit_system(),
                              m_distributions.size(), m_packing_dimensions); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMFEvoR::write_variables(const std::set<std::string> &vars, const PIO& nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "enhancement_factor")) {
    ierr = m_enhancement_factor.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "distributions")) {
    ierr = save_distributions(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "recrystallizations")) {
    ierr = save_diagnostics(nc); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMFEvoR::allocate() {
  PetscErrorCode ierr;

  // SIAFD diffusive flux computation requires stencil width of 2 (we
  // traverse the grid *with ghosts* and we need staggered grid offset
  // data at each point)
  const unsigned int stencil_width = 2;

  ierr = m_enhancement_factor.create(grid, "enhancement_factor", WITH_GHOSTS,
                                     stencil_width); CHKERRQ(ierr);
  ierr = m_enhancement_factor.set_attrs("diagnostic", // i.e. not needed to re-start the model
                                        "flow law enhancement factor",
                                        "1", // dimensionless
                                        "" /* no standard name */); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMFEvoR::init(Vars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
                    "* Initializing the Fabric Evolution with Recrystallization"
                    " (FEvoR) model...\n"); CHKERRQ(ierr);

  // set enhancement factor to something reasonable
  ierr = m_enhancement_factor.set(config.get("sia_enhancement_factor")); CHKERRQ(ierr);
  // make enhancement factor available to other PISM components:
  ierr = vars.add(m_enhancement_factor); CHKERRQ(ierr); //IceModelVec

  m_enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (m_enthalpy == NULL) {
    SETERRQ(grid.com, 1, "enthalpy field is not available");
  }

  // FIXME: load packing dimensions here

  // It would be nice to be able to allocate these in
  // PISMFEvoR::allocate() or in the constructor, but pism::Vars is
  // not available there...
  if (m_pressure == NULL) {
    m_pressure = new PSB_pressure(m_stress_balance, grid, vars);
    assert(m_pressure != NULL);
  }

  if (m_tauxz == NULL) {
    m_tauxz = new PSB_tauxz(m_stress_balance, grid, vars);
    assert(m_tauxz != NULL);
  }

  if (m_tauyz == NULL) {
    m_tauyz = new PSB_tauyz(m_stress_balance, grid, vars);
    assert(m_tauyz != NULL);
  }

  bool i_set = false;
  std::string input_file;
  ierr = OptionsString("-i", "input file name", input_file, i_set); CHKERRQ(ierr);

  if (i_set) {
    ierr = load_distributions(input_file); CHKERRQ(ierr);
  } else {
    ierr = set_initial_distribution_parameters(); CHKERRQ(ierr);
  }


  return 0;
}


/** Create an initial state (cloud of "particles").
 *
 * Use this for testing or bootstrapping only.
 *
 * @return 0 on success
 */
PetscErrorCode PISMFEvoR::set_initial_distribution_parameters() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "  Setting initial distribution parameters...\n"); CHKERRQ(ierr);


  unsigned int n_particles = (grid.Mz-2)*(grid.Mx-1);
  
  // Initialize distributions
  assert(m_packing_dimensions.size() == 3);
  double w_i = -3.0; // This makes a weak bi-polar (single maximum)
  FEvoR::Distribution d_i(m_packing_dimensions, w_i);

  m_distributions.resize(n_particles, FEvoR::Distribution(m_packing_dimensions, w_i));

  // Initialize particle positions and corresponding enhancement factors
  m_p_x.resize(n_particles);
  m_p_y.resize(n_particles);
  m_p_z.resize(n_particles);
  m_p_e.resize(n_particles);
  
  for (unsigned int zz = 0; zz < grid.Mz-2; ++zz) {
    for (unsigned int xx = 0; xx < grid.Mx-1; ++xx) {
      unsigned int i = xx + zz*(grid.Mx-1);
      
      m_p_x[i] = grid.x[xx];
      m_p_y[i] = 0.0;
      m_p_z[i] = (grid.zlevels[zz]+grid.zlevels[zz+1])/2.0;
      m_p_e[i] = 1.0;
    }
  }

  return 0;
}

/** Load distributions and particle positions from a file specified using the "-i" option.
 *
 * Uses the last time record.
 *
 * @return 0 on success
 */
PetscErrorCode PISMFEvoR::load_distributions(const std::string &input_file) {

  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "  Reading distribution data from %s...\n",
                    input_file.c_str()); CHKERRQ(ierr);

  PIO nc(grid, "guess_mode");
  ierr = nc.open(input_file, PISM_READONLY); CHKERRQ(ierr);
  {
    unsigned int n_records = 0;
    ierr = nc.inq_nrecords(n_records); CHKERRQ(ierr);
    unsigned int last_record = n_records - 1;

    // Get packing dimensions:
    std::vector<double> pd;
    ierr = nc.get_att_double("PISM_GLOBAL", "packing_dimensions", pd); CHKERRQ(ierr);

    if (pd.size() != 3) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: cannot load FEvoR distributions from %s.\n",
                         input_file.c_str()); CHKERRQ(ierr);
      PISMEnd();
    }

    for (unsigned int k = 0; k < 3; ++k) {
      m_packing_dimensions[k] = pd[k];
    }

    // Get the number of distributions:
    unsigned int n_particles = (m_packing_dimensions[0] *
                                m_packing_dimensions[1] *
                                m_packing_dimensions[2]);

    // Create the isotropic distribution:
    m_d_iso = FEvoR::Distribution(m_packing_dimensions, 0.0);

    // Fill all distributions with copies of the isotropic one:
    m_distributions = std::vector<FEvoR::Distribution>(n_particles, m_d_iso);

    for (unsigned int k = 0; k < n_particles; ++k) {
      ierr = fevor_load_distribution(nc, k, last_record, m_distributions[k]); CHKERRQ(ierr);
    }

    // Load particle positions:
    m_p_x.resize(n_particles, 0);
    m_p_y.resize(n_particles, 0);
    m_p_z.resize(n_particles, 0);
    ierr = fevor_load_particle_positions(nc, last_record, m_p_x, m_p_y, m_p_z); CHKERRQ(ierr);

    // Resize storage for enhancement factors:
    m_p_e.resize(n_particles, 0);
  }
  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

/** Save distributions and particle positions to a file
 *
 * Uses the last time record.
 *
 * @param nc file to save to
 *
 * @return 0 on success
 */
PetscErrorCode PISMFEvoR::save_distributions(const PIO &nc) {

  PetscErrorCode ierr;

  unsigned int n_records = 0;
  ierr = nc.inq_nrecords(n_records); CHKERRQ(ierr);
  unsigned int last_record = n_records - 1;

  unsigned int n_particles = m_distributions.size();

  for (unsigned int k = 0; k < n_particles; ++k) {
    ierr = fevor_save_distribution(nc, k, last_record, m_distributions[k]); CHKERRQ(ierr);
  }

  // Save particle positions:
  ierr = fevor_save_particle_positions(nc, last_record, m_p_x, m_p_y, m_p_z); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMFEvoR::save_diagnostics(const PIO &nc) {
  PetscErrorCode ierr;

  unsigned int n_particles = m_distributions.size();

  std::vector<double> tmp_m(n_particles), tmp_p(n_particles);
  for (unsigned int k = 0; k < n_particles; ++k) {
    tmp_m[k] = m_n_migration_recrystallizations[k];
    tmp_p[k] = m_n_polygonization_recrystallizations[k];
  }

  unsigned int n_records = 0;
  ierr = nc.inq_nrecords(n_records); CHKERRQ(ierr);
  unsigned int last_record = n_records - 1;

  ierr = fevor_save_recrystallization_numbers(nc, last_record, tmp_m, tmp_p); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
