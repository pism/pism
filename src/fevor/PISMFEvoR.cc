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
#include <iomanip>
#include <math.h>

#include <fevor_distribution.hh>
#include <vector_tensor_operations.hh>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include "PISMFEvoR.hh"
#include "base/util/PISMConfig.hh"
#include "base/util/PISMVars.hh"
#include "base/stressbalance/PISMStressBalance_diagnostics.hh"
#include "base/enthalpyConverter.hh"
#include "base/util/error_handling.hh"

#include "FEvoR_IO.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMTime.hh"
#include "base/util/io/PIO.hh"
#include "base/util/IceModelVec.hh"
#include "base/util/MaxTimestep.hh"

#include "base/util/pism_memory.hh"
using PISM_SHARED_PTR_NSPACE::dynamic_pointer_cast;

#include<iostream>

namespace pism {

  PISMFEvoR::PISMFEvoR(IceGrid::ConstPtr g, stressbalance::StressBalance *sb)
  : Component_TS(g),
    m_stress_balance(sb),
      m_packing_dimensions(std::vector<unsigned int>(3, 16)),
      m_d_iso(m_packing_dimensions, 0.0)
  {
 
  unsigned int n_particles = (unsigned int)m_config->get_double("fevor_n_particles");
                
  m_distributions = std::vector<FEvoR::Distribution>(n_particles, m_d_iso);
  m_p_avg_stress = std::vector<double>(n_particles * 9, 0.0);
  m_p_avg_temp = std::vector<double>(n_particles, 0.0);

  m_EC = m_grid->ctx()->enthalpy_converter();
  assert(m_stress_balance != NULL);

  //  allocate(); 

  const unsigned int stencil_width = 2;

  m_enhancement_factor = new IceModelVec3;
  m_enhancement_factor->create(m_grid, "enhancement_factor", WITH_GHOSTS,
                                     stencil_width);
  m_enhancement_factor->set_attrs("diagnostic", // i.e. not needed to re-start the model
                                        "flow law enhancement factor",
                                        "1", // dimensionless
                                        "" /* no standard name */);
  last_update_time = m_grid->ctx()->time()->start(); // Actually this should be overwritten by init with time of last update when loading from restart file. // FLO // TODO // ERROR // FIXME

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

MaxTimestep PISMFEvoR::max_timestep_impl(double t) {
  // uses current calander definition of a year
  //m_grid->ctx()->time()->current()
  const double years_per_second = m_grid->ctx()->time()->convert_time_interval(1., "years");
  
  const double m_start_time = m_grid->ctx()->time()->start();
  
  double fevor_step = m_config->get_double("fevor_step");
  fevor_step = fevor_step / years_per_second;

  double fevor_dt = (fevor_step - (t - last_update_time));
  assert (fevor_dt > 0 ) ; 
  const double epsilon = 1.0; // 1 second tolerance
  if (fevor_dt > epsilon) { 
    MaxTimestep dt = fevor_dt;
    return dt;
  }  else
    {
      MaxTimestep dt; 
      return dt;
    }
}

void PISMFEvoR::update_impl(double t, double dt) {


  m_t = t;
  m_dt = dt;

  // determine whether or not to run FEvoR
  // -------------------------------------
  bool step_flag = false;
  const double epsilon = 1.0; // 1 second tolerance
  // uses current calander definition of a year
  double years_per_second = m_grid->ctx()->time()->convert_time_interval(1., "years");
  
  double fevor_step = (double)m_config->get_double("fevor_step");
  fevor_step = fevor_step / years_per_second;
  double fevor_dt = m_t + m_dt - last_update_time;
  
  if (fevor_dt >= fevor_step*.9+epsilon) {
    step_flag = true;
    last_update_time= m_t+m_dt;
  }
  if (step_flag)  m_log->message(4,
                    "\n actually doing update_impl in FEvoR()\n"); 
  
  const IceModelVec3
    &u = m_stress_balance->velocity_u(),
    &v = m_stress_balance->velocity_v(),
    &w = m_stress_balance->velocity_w();

  assert(m_pressure != NULL);
  IceModelVec::Ptr pressure  = m_pressure->compute();
  IceModelVec3::Ptr pressure3 = dynamic_pointer_cast<IceModelVec3,IceModelVec>(pressure);


  assert(m_tauxz != NULL);
  IceModelVec::Ptr   tauxz = m_tauxz->compute(); 
  IceModelVec3::Ptr tauxz3    = dynamic_pointer_cast<IceModelVec3,IceModelVec>(tauxz);

  assert(m_tauyz != NULL);
  IceModelVec::Ptr tauyz = m_tauyz->compute();
  IceModelVec3::Ptr tauyz3 = dynamic_pointer_cast<IceModelVec3,IceModelVec>(tauyz);

  {
    IceModelVec::AccessList list;
    list.add(u);
    list.add(v);
    list.add(w);
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
      assert(0.0 <= m_p_z[i] && m_p_z[i] <= m_grid->Lz());
      assert(m_grid->x().front() <= m_p_x[i] && m_p_x[i] <= m_grid->x().back());
      assert(m_grid->y().front() <= m_p_y[i] && m_p_y[i] <= m_grid->y().back());

      evaluate_at_point(*pressure3, m_p_x[i], m_p_y[i], m_p_z[i], P);    
      evaluate_at_point(*tauxz3,    m_p_x[i], m_p_y[i], m_p_z[i], txz);  
      evaluate_at_point(*tauyz3,    m_p_x[i], m_p_y[i], m_p_z[i], tyz);  

      // Make sure that the pressure is positive. (Zero pressure
      // probably means that the current particle exited the ice blob.
      // This is a failure of the PISMFEvoR logic: such particles need
      // to be discarded and replaced by new ones *within* the ice.)
      assert(P >= 0.0);

      evaluate_at_point(*m_enthalpy, m_p_x[i], m_p_y[i], m_p_z[i], E); 
      T = m_EC->temperature(E, P); 
     
      //m_p_avg_temp[i] += T*m_dt;
      m_p_avg_temp[i] = T;
      /* Indexing: {0, 1, 2,
       *            3, 4, 5,
       *            6, 7, 8}
       */
      //m_p_avg_stress[i*9 + 0] = m_p_avg_stress[i*9 + 4] = m_p_avg_stress[i*9 + 8] += -P*m_dt;
      //m_p_avg_stress[i*9 + 2] = m_p_avg_stress[i*9 + 6] += txz*m_dt; 
      //m_p_avg_stress[i*9 + 5] = m_p_avg_stress[i*9 + 7] += tyz*m_dt; 
      m_p_avg_stress[i*9 + 0] = m_p_avg_stress[i*9 + 4] = m_p_avg_stress[i*9 + 8] = -P;
      m_p_avg_stress[i*9 + 2] = m_p_avg_stress[i*9 + 6] = txz; 
      m_p_avg_stress[i*9 + 5] = m_p_avg_stress[i*9 + 7] = tyz; 



      if (step_flag) {
        double fevor_begin = m_t + m_dt - fevor_step;
        //double temp = m_p_avg_temp[i]/fevor_step;
        double temp = m_p_avg_temp[i];
        m_p_avg_temp[i] = 0.;
        /* Indexing: {0, 1, 2,
         *            3, 4, 5,
         *            6, 7, 8}
         */
        std::vector<double> stress(9,0.0);
        for (unsigned int s = 0; s < 9; ++s){
          //stress[s] = m_p_avg_stress[i*9+s]/fevor_step;
          stress[s] = m_p_avg_stress[i*9+s];
          m_p_avg_stress[i*9+s] = 0.;
        }

        std::vector<double> bulkEdot(9, 0.0);
        std::vector<double> bulkM(81, 0.0);
        
        assert(temp != 0.);
        // http://en.wikipedia.org/wiki/Step_in_Time
        bulkM = m_distributions[i].stepInTime(temp, stress, last_update_time, fevor_dt,
                                              m_n_migration_recrystallizations[i],
                                              m_n_polygonization_recrystallizations[i],
                                              bulkEdot);
        
        std::vector<double> bulkEdot_iso(9, 0.0);
        std::vector<double> bulkM_iso(81, 0.0);
        m_d_iso = FEvoR::Distribution(m_packing_dimensions, 0.0); // make sure this is a clean distribution. 
        bulkM_iso = m_d_iso.stepInTime(temp, stress, 0., 0., bulkEdot_iso);
        if (bulkEdot_iso[2] != 0.0) {
          m_p_e[i] = std::abs(bulkEdot[2] / bulkEdot_iso[2]);
        } else {
          // If a particle is outside the ice blob, then pressure == 0
          // and bulkEdot_iso == 0, so we need this special case.
          m_p_e[i] = 1.0;
        }
        
        // // some bounds for the enhancement factor -- in shear only
        // if (m_p_e[i] < 1.0) {
        //   m_p_e[i] = 1.0;
        // } else if (m_p_e[i] > 10.0) {
        //   m_p_e[i] = 10.0;
        //   // upper bound.
        // }
      }
      
    } // end of the for-loop over particles
    
    if (step_flag) {
      // set the enhancement factor for every grid point from our particle cloud
      pointcloud_to_grid(m_p_x, m_p_z, m_p_e, *m_enhancement_factor); 

      m_enhancement_factor->update_ghosts(); 
    }

    // update particle positions -- don't want to update until gridded
    // enhancement factor is updated (if it is)
    for (unsigned int i = 0; i < n_particles; ++i) {
      double p_u = 0.0,
             p_v = 0.0,
             p_w = 0.0;

      evaluate_at_point(u, m_p_x[i], m_p_y[i], m_p_z[i], p_u);   
      evaluate_at_point(v, m_p_x[i], m_p_y[i], m_p_z[i], p_v);   
      evaluate_at_point(w, m_p_x[i], m_p_y[i], m_p_z[i], p_w);   

      update_particle_position(m_p_x[i], m_p_y[i], m_p_z[i], p_u, p_v, p_w, m_dt); 
    }
  }

  //  delete & (* pressure);
  //  delete & (* tauxz);
  //  delete & (* tauyz);

    m_log->message(4,
                    "\n END update_impl in FEvoR()\n"); 

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
void PISMFEvoR::update_particle_position(double &x, double &y, double &z,
                                                   double u, double v, double w,
                                                   double dt) {

  // basic Euler method. Assuming u, v, w are 'good'
  x += u*dt;
  y += v*dt; // probably not needed and v should be zero in 2D flow line model
  z += w*dt; // w should be zero

  if (z > m_grid->Lz()) {
    z = m_grid->Lz(); 
  } else if (z < 0.0) { // WE ARE IN THE BEDROCK! DO SOMETHING !  FLO ! TODO
    z = 0.0;
  }
   
  // assuming periodic x... 
  if (x < m_grid->x().front()) {
    x += m_grid->x().back()-m_grid->x().front();
  } else if (x > m_grid->x().back()) {
    x -= m_grid->x().back()-m_grid->x().front();
  }
  
  // assuming periodic y... 
  if (y < m_grid->y().front()) {
    y += m_grid->y().back()-m_grid->y().front();
  } else if (y > m_grid->y().back()) {
    y -= m_grid->y().back()-m_grid->y().front();
  }

  // check if the point (x,y,z) is within the domain:
  assert(0.0 <= z && z <= m_grid->Lz());
  // PISM's horizontal grid is cell-centered, so (-m_grid->Lx(),
  // -m_grid->Ly) is actually outside of the convex hull of all grid
  // points.
  assert(m_grid->x().front() <= x && x <= m_grid->x().back());
  assert(m_grid->y().front() <= y && y <= m_grid->y().back());

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
void PISMFEvoR::evaluate_at_point(const IceModelVec3 &input,
                                            double x, double y, double z,
                                            double &result) {

  // compute indexes for accessing neighboring columns:
  int I_left = 0, I_right = 0, J_bottom = 0, J_top = 0;
  m_grid->compute_point_neighbors(x, y, I_left, I_right, J_bottom, J_top);

  // compute index of the vertical level just below z:
  unsigned int K = 0;
  // K + 1 (used below) should be at most Mz() - 1
  while (K + 1 < m_grid->Mz() - 1 && m_grid->z(K + 1) < z) {
    K++;
  }

  // get pointers to neighboring columns:
  const double* column[4] = {NULL, NULL, NULL, NULL};
  column[0] = input.get_column(I_left,  J_bottom);
  column[1] = input.get_column(I_right, J_bottom);
  column[2] = input.get_column(I_right, J_top);
  column[3] = input.get_column(I_left,  J_top);

  // compute interpolation weights
  std::vector<double> weights = m_grid->compute_interp_weights(x, y);
  assert(K < m_grid->Mz());
  const double z_weight = (z - m_grid->z(K)) / (m_grid->z(K+1) - m_grid->z(K));

  // interpolate
  result = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    // vertical interpolation:
    double f_i = column[i][K] + z_weight * (column[i][K+1] - column[i][K]);
    // horizontal interpolation:
    result += weights[i] * f_i;
  }

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
void PISMFEvoR::pointcloud_to_grid(const std::vector<double> &x,
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

  //FLO // WEIGHTS FOR INTERPOLATION
  const double
    wx=1./m_grid->dx(), 
    wz=1./m_grid->z(1);
  
  // get the points for our convex hull
  for (unsigned int pn = 0; pn < n_particles; ++pn) {
    Point p(x[pn]*wx,z[pn]*wz);
    D_TRI.insert( p );
    function_values.insert( std::make_pair(p, values[pn]) );
  }

  IceModelVec::AccessList list(result);

  for (int i=m_grid->xs()+1; i<m_grid->xs()+m_grid->xm()-1; ++i) {
    for (unsigned int k=1; k<m_grid->Mz(); ++k) {
      Point INTERP( m_grid->x(i)*wx,m_grid->z(k)*wz );

      // make a vector of the iterpolation point and type
      Point_coordinate_vector coord;
      The_neighbors norm = CGAL::natural_neighbor_coordinates_2 (D_TRI, INTERP, std::back_inserter(coord) );
      
      Field_type res;
      // make sure point is in convex hull!
      if(!norm.third) {
        // FIXME this is dumb logic. Works for only the Enhancement factor now!!
        if (m_grid->z(k) > m_grid->z(m_grid->Mz()-1)/2.0) {
          // above the middle
          res = Field_type(1.0); // isotropic.
        } else {
          res = Field_type(2.0);// Field_type( *std::max_element(values.begin(), values.end()) ); // max.
	  m_log->message(2,
			 "One point outside of convex hull for enthalpy.\n"
			 "Setting E to 2.0.\n"
			 "x-value of missing point %f\n"
			 "z-value of missing point %f", m_grid->x(k) ,m_grid->z(k));
	} 
      } else {
        res =  CGAL::linear_interpolation (coord.begin(), coord.end(), norm.second, Value_access(function_values));
      }
      
      // set result for all y grid points at INTERP(x,z)
      for (int j=m_grid->ys(); j<m_grid->ys()+m_grid->ym(); ++j) {
        result(i,j,k) = double(res);
      }
    }
      for (int j=m_grid->ys(); j<m_grid->ys()+m_grid->ym(); ++j) {
        result(i,j,0) = result(i,j,1);
      }
  }
  
  // deal with side boundaries. 
  for (int j=m_grid->ys(); j<m_grid->ys()+m_grid->ym(); ++j) {
    for (unsigned int k=0; k<m_grid->Mz(); ++k) {
      result(m_grid->xs(),j,k) = result(m_grid->xs()+1,j,k);
      result(m_grid->xs()+m_grid->xm()-1,j,k) = result(m_grid->xs()+m_grid->xm()-2,j,k);
    }
  } 
  
}

void PISMFEvoR::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  if (keyword != "none") {
    result.insert(m_enhancement_factor->metadata().get_string("short_name"));
    result.insert("distributions");
  }
}

void PISMFEvoR::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                           IO_Type nctype) {
  if (set_contains(vars, "enhancement_factor")) {
     m_enhancement_factor->define(nc, nctype);
  }

  if (set_contains(vars, "distributions") or set_contains(vars, "recrystallizations")) {
     fevor_prepare_file(nc, m_grid->ctx()->unit_system(),
                              m_distributions.size(), m_packing_dimensions);
  }

  
}

void PISMFEvoR::write_variables_impl(const std::set<std::string> &vars, const PIO& nc) {
  if (set_contains(vars, "enhancement_factor")) {
     m_enhancement_factor->write(nc);
  }

  if (set_contains(vars, "distributions")) {
     save_distributions(nc);
  }

  if (set_contains(vars, "recrystallizations")) {
     save_diagnostics(nc);
  }

  
}

void PISMFEvoR::allocate() {

  // SIAFD diffusive flux computation requires stencil width of 2 (we
  // traverse the grid *with ghosts* and we need staggered grid offset
  // data at each point)

  
}

void PISMFEvoR::init() {
   m_log->message(2,
                    "* Initializing the Fabric Evolution with Recrystallization"
                    " (FEvoR) model...\n");

  // set enhancement factor to something reasonable
  m_enhancement_factor->set(m_config->get_double("sia_enhancement_factor")); 
  // make enhancement factor available to other PISM components:
  //   m_grid->variables().add(m_enhancement_factor);  //IceModelVec // FLO ! // NEEDS TO HAPPEN // ERROR // TODO

  m_enthalpy =  m_grid->variables().get_3d_scalar("enthalpy");
  if (m_enthalpy == NULL) {
    throw RuntimeError("enthalpy field is not available");
  }

  m_log->message(2,
                    "sorting stressbalance stuff"
                    " (FEvoR) model...\n");

  // FIXME: load packing dimensions here

  // It would be nice to be able to allocate these in
  // PISMFEvoR::allocate() or in the constructor, but pism::Vars is
  // not available there...
  if (m_pressure == NULL) {
    m_pressure = new stressbalance::PSB_pressure(m_stress_balance);
    assert(m_pressure != NULL);
  }

  if (m_tauxz == NULL) {
    m_tauxz = new stressbalance::PSB_tauxz(m_stress_balance);
    assert(m_tauxz != NULL);
  }

  if (m_tauyz == NULL) {
    m_tauyz = new stressbalance::PSB_tauyz(m_stress_balance);
    assert(m_tauyz != NULL);
  }
   m_log->message(2,
                    "Going to load distributions"
                    " (FEvoR) model...\n");

  std::string input_file;
  int start = 0;
  bool boot = false;
  bool use_input_file = find_pism_input(input_file, boot, start);

  if (use_input_file && ! boot) {
     load_distributions(input_file);
  } else {
     set_initial_distribution_parameters();
  }


   m_log->message(2,
                    "Loaded distributions"
                    " (FEvoR) model...\n");
  
}
  IceModelVec3 * PISMFEvoR::get_enhancement_factor(){
    return (m_enhancement_factor);
  }

/** Create an initial state (cloud of "particles").
 *
 * Use this for testing or bootstrapping only.
 *
 * @return 0 on success
 */
void PISMFEvoR::set_initial_distribution_parameters() {

   verbPrintf(2, m_grid->com,
                    "  Setting initial distribution parameters...\n");
  
  unsigned int n_particles = (unsigned int)m_config->get_double("fevor_n_particles");
  assert( n_particles == (m_grid->Mz()-1)*(m_grid->Mx()-1) );
  
  // Initialize distributions
  assert(m_packing_dimensions.size() == 3);
  double w_i = -3.0; // This makes a weak bi-polar (single maximum)
  FEvoR::Distribution d_i(m_packing_dimensions, w_i);

  m_distributions.resize(n_particles, d_i);

  // Initialize particle positions and corresponding enhancement factors
  m_p_x.resize(n_particles);
  m_p_y.resize(n_particles);
  m_p_z.resize(n_particles);
  m_p_e.resize(n_particles);
  
  for (unsigned int zz = 0; zz < m_grid->Mz()-1; ++zz) {
    for (unsigned int xx = 0; xx < (unsigned int)m_grid->Mx()-1; ++xx) {
      unsigned int i = xx + zz*(m_grid->Mx()-1);
      unsigned int invZ = m_grid->Mz()-zz-1;
      
      m_distributions[i].generateWatsonAxes(w_i*double(invZ)); 
        // increase fabric strength with depth
      
      m_p_x[i] = m_grid->x(xx);
      m_p_y[i] = 0.0;
      m_p_z[i] = (m_grid->z(zz)+m_grid->z(zz+1))/2.0;
      m_p_e[i] = 1.0;
    }
  }

  
}

/** Load distributions and particle positions from a file specified using the "-i" option.
 *
 * Uses the last time record.
 *
 * @return 0 on success
 */
void PISMFEvoR::load_distributions(const std::string &input_file) {


   verbPrintf(2, m_grid->com,
                    "  Reading distribution data from %s...\n",
                    input_file.c_str());

  PIO nc(m_grid->com, "guess_mode");
   nc.open(input_file, PISM_READONLY);
  {
    unsigned int n_records = nc.inq_nrecords();
    unsigned int last_record = n_records - 1;

    // Get packing dimensions:
    std::vector<double> pd = nc.get_att_double("PISM_GLOBAL", "packing_dimensions");

    if (pd.size() != 3) {
      throw RuntimeError::formatted( "PISM ERROR: cannot load FEvoR distributions from %s.\n",
                         input_file.c_str());
    }

    for (unsigned int k = 0; k < 3; ++k) {
      m_packing_dimensions[k] = pd[k];
    }

    // Get the number of distributions:
    unsigned int n_particles = 0;
    n_particles = nc.inq_dimlen("distribution_index");

    // Create the isotropic distribution:
    m_d_iso = FEvoR::Distribution(m_packing_dimensions, 0.0);

    // Fill all distributions with copies of the isotropic one:
    m_distributions = std::vector<FEvoR::Distribution>(n_particles, m_d_iso);

    for (unsigned int k = 0; k < n_particles; ++k) {
       fevor_load_distribution(nc, k, last_record, m_distributions[k]);
    }

    // Load particle positions:
    m_p_x.resize(n_particles, 0);
    m_p_y.resize(n_particles, 0);
    m_p_z.resize(n_particles, 0);
     fevor_load_particle_positions(nc, last_record, m_p_x, m_p_y, m_p_z);

    // Resize storage for enhancement factors:
    m_p_e.resize(n_particles, 0);
    // Resize storage for average stress and temperature between fevor steps
    m_p_avg_stress.resize(n_particles * 9, 0.0);
    m_p_avg_temp.resize(n_particles, 0.0);

  }
   nc.close();

  
}

/** Save distributions and particle positions to a file
 *
 * Uses the last time record.
 *
 * @param nc file to save to
 *
 * @return 0 on success
 */
void PISMFEvoR::save_distributions(const PIO &nc) {


  unsigned int n_records = 0;
  n_records = nc.inq_nrecords();
  unsigned int last_record = n_records - 1;

  unsigned int n_particles = m_distributions.size();

  for (unsigned int k = 0; k < n_particles; ++k) {
     fevor_save_distribution(nc, k, last_record, m_distributions[k]);
  }

  // Save particle positions:
  fevor_save_particle_positions(nc, last_record, m_p_x, m_p_y, m_p_z, m_p_e);

  
}

void PISMFEvoR::save_diagnostics(const PIO &nc) {

  unsigned int n_particles = m_distributions.size();

  std::vector<double> tmp_m(n_particles), tmp_p(n_particles);
  for (unsigned int k = 0; k < n_particles; ++k) {
    tmp_m[k] = m_n_migration_recrystallizations[k];
    tmp_p[k] = m_n_polygonization_recrystallizations[k];
  }

  unsigned int n_records = 0;
  n_records = nc.inq_nrecords();
  unsigned int last_record = n_records - 1;

   fevor_save_recrystallization_numbers(nc, last_record, tmp_m, tmp_p);

  
}

} // end of namespace pism
