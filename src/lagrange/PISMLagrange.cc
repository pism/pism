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
#include <list>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <random>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include <petscdmda.h>
#include "petscdmda.h"   


#include "PISMLagrange.hh"
#include "base/util/PISMConfig.hh"
#include "base/util/PISMVars.hh"
#include "base/stressbalance/PISMStressBalance_diagnostics.hh"
#include "base/enthalpyConverter.hh"

#include "lagrange_IO.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMTime.hh"
#include "base/util/Profiling.hh"
#include "base/util/io/PIO.hh"

#include<iostream>
#include <unistd.h>
#include <mpi.h>

namespace pism {

  PISMLagrange::PISMLagrange(IceGrid::ConstPtr g, stressbalance::StressBalance *stress_balance)
  : Component_TS(g), m_stress_balance(stress_balance)
{
 
  assert(m_stress_balance != NULL);

  allocate();
}

PISMLagrange::~PISMLagrange() {
}

  /** Make timesteps hit seed times
   * 
   * Basically there's no real limitation on timestep size. 
   * 2D CFL should ensure reasonable timesteps, so we don't 
   * have to care about that.
   *
   */
  MaxTimestep PISMLagrange::max_timestep_impl(double t) {
  double dt = 0 ;
  if (seed_times.size() == 0 || t >= seed_times.back()) {
    return MaxTimestep();
  }

  if (t < *next_seed) 
  dt = *next_seed - t;

  // now make sure that we don't end up taking a time-step of less than 1
  // second long
  if (dt < 1.0) {
    if (next_seed != seed_times.end()-1 ) {
      return MaxTimestep(*(next_seed+1) - t);
    } else {
      return MaxTimestep();
    }
  } else {
    return MaxTimestep(dt);
  }
}

  
void PISMLagrange::update_impl(double t, double dt) {
  m_log->message(4,
		   "Starting tracer update\n");
  m_t = t;
  m_dt = dt;

  m_log->message(5,
		   "Removing flying tracers\n");
  remove_flying_tracers();

  m_log->message(5,
		   "(possibly) seeding\n");
  if (check_seed())
    seed();

  m_log->message(5,
		   "getting 3d fields\n");
  const IceModelVec3
    &u = m_stress_balance->velocity_u(),
    &v = m_stress_balance->velocity_v(),
    &w = m_stress_balance->velocity_w();
  {
    IceModelVec::AccessList list;
    list.add(u);
    list.add(v);
    list.add(w);

    m_log->message(5,
		   "going through the list\n");
    for (std::list<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
      const double
	p_u = evaluate_at_point(u, it->x, it->y, it->z),
	p_v = evaluate_at_point(v, it->x, it->y, it->z),
	p_w = evaluate_at_point(w, it->x, it->y, it->z);
      update_particle_position(it->x, it->y, it->z, p_u, p_v, p_w, m_dt);
    }
  }

  m_log->message(5,
		   "removing bedrock tracers\n");
  remove_bedrock_tracers();

  m_log->message(5,
		   "shipping tracers\n");
  ship_tracers(); // Needs to be done before positions get messed up by periodic boundary alignment

  m_log->message(5,
		   "making things periodic\n");
  /// @todo: Check if boundaries are supposed to be periodic before pushing particles around // TODO // 
  for (std::list<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {

    // assuming periodic x... 
    if (it->x < m_grid->x().front() - m_grid->dx()/2.) {
      it->x += m_grid->Lx()*2.;
    } else if (it->x > m_grid->x().back() + m_grid->dx()/2.) {
      it->x -= m_grid->Lx()*2.;
    }
    
    // assuming periodic y... 
    if (it->y < m_grid->y().front() - m_grid->dy()/2.) {
      it->y += m_grid->Ly()*2.;
    } else if (it->y > m_grid->y().back() + m_grid->dy()/2.) {
      it->y -= m_grid->Ly()*2.;
    }
  }
  m_log->message(4,
		   "done updating tracers\n");
}

/**
 * Evaluate new particle position from velocity
 *
 * @param x, y, z coordinates of a point within the domain
 * @param u, v, w velocities of x, y, z
 * @param dt time step length
 *
 */
void PISMLagrange::update_particle_position(double &x, double &y, double &z,
                                                   double u, double v, double w,
                                                   double dt) {

  // basic Euler method. Assuming u, v, w are 'good'
  x += u*dt;
  y += v*dt; // probably not needed and v should be zero in 2D flow line model
  z += w*dt; // w should be zero

}

/**
 * Evaluate a 3D field at a given point.
 *
 * @param input a 3D field
 * @param x, y, z coordinates of a point within the domain
 *
 */
double PISMLagrange::evaluate_at_point(const IceModelVec3 &input,
                                            const double x, const double y, const double z) {
  //PetscErrorCode ierr;

  // compute indexes for accessing neighboring columns:
  int I_left = 0, I_right = 0, J_bottom = 0, J_top = 0;
  m_grid->compute_point_neighbors(x, y, I_left, I_right, J_bottom, J_top);


  
  // compute index of the vertical level just below z:
  unsigned int K = 0;
  // K + 1 (used below) should be at most Mz - 1
  while (K + 1 < m_grid->Mz() - 1 && m_grid->z(K + 1) < z) {
    K++;
   }


  unsigned int index_i[4] = {I_left, I_right, I_right, I_left};
  unsigned int index_j[4] = {J_bottom, J_bottom, J_top, J_top};


  // compute interpolation weights
  std::vector<double> weights = m_grid->compute_interp_weights(x, y);
  assert(K+1 < m_grid->Mz());
  const double z_weight = (z - m_grid->z(K)) / (m_grid->z(K+1) - m_grid->z(K));

  // interpolate
  double result = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    // vertical interpolation:
    const double below = input(index_i[i], index_j[i], K),
                 above = input(index_i[i], index_j[i], K+1);
    double f_i =  below + z_weight * (above - below);
    // horizontal interpolation:
    result += weights[i] * f_i;
  }
  return result;
}




/**
 * Evaluate a 2D field at a given point.
 *
 * @param input a 2D field
 * @param x, y coordinates of a point within the domain
 *
 */
double PISMLagrange::evaluate_at_point(const IceModelVec2S &input,
                                            double x, double y) {
  // compute indexes for accessing neighboring columns:
  int I_left = 0, I_right = 0, J_bottom = 0, J_top = 0;
  m_grid->compute_point_neighbors(x, y, I_left, I_right, J_bottom, J_top);


  unsigned int index_i[4] = {I_left, I_right, I_right, I_left};
  unsigned int index_j[4] = {J_bottom, J_bottom, J_top, J_top};


  // compute interpolation weights
  std::vector<double> weights = m_grid->compute_interp_weights(x, y);
  // interpolate
  double result = 0.0;
  for (unsigned int i = 0; i < 4; ++i) {
    // vertical interpolation:
    const double contribution = input(index_i[i], index_j[i]);
    result += weights[i] * contribution;
  }
  return result;
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
 */
void PISMLagrange::pointcloud_to_grid(const std::vector<double> &x,
						const std::vector<double> &y,
                                             const std::vector<double> &z,
                                             const std::vector<double> &values,
                                             IceModelVec3 &result) {

  assert(x.size() ==y.size() && x.size() ==z.size() && z.size() == values.size());

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

  for (int i=m_grid->xs(); i<m_grid->xs()+m_grid->xm(); ++i) {
    for (unsigned int k=0; k<m_grid->Mz(); ++k) {
      Point INTERP( m_grid->x(i),m_grid->z(k) );

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
          res = Field_type( *std::max_element(values.begin(), values.end()) ); // max.
        } 
      } else {
        res =  CGAL::linear_interpolation (coord.begin(), coord.end(), norm.second, Value_access(function_values));
      }
      
      // set result for all y grid points at INTERP(x,z)
      for (int j=m_grid->ys(); j<m_grid->ys()+m_grid->ym(); ++j) {
        result(i,j,k) = double(res);
      }
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

void PISMLagrange::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
}

void PISMLagrange::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                           IO_Type nctype) {
    int ps = particles.size(), as = 0;
    MPI_Allreduce(&ps, &as, 1, MPI_INT, MPI_SUM, m_grid->com);
    if (as == 0 ) as = 1 ; // Otherwise defining the variables will fail.
    lagrange_prepare_file(nc, as);
}

  void PISMLagrange::write_variables_impl(const std::set<std::string> &vars,const  PIO& nc) {

    save_particle_positions(nc);
 
}

void  PISMLagrange::allocate() {


}


void PISMLagrange::init() {
  //Define the particle type -- probably only needs to be done once.  
  MPI_Datatype oldtypes[2]; 
  int          blockcounts[2];
  MPI_Aint    offsets[2], extent, tester;
  offsets[0] = 0;
  oldtypes[0] = MPI_DOUBLE;
  blockcounts[0] = 3;
  /* Setup description of the 2 MPI_INT fields n, type */
  /* Need to first figure offset by getting size of MPI_FLOAT */
  MPI_Type_get_extent(MPI_DOUBLE, &tester, &extent);
  offsets[1] = 3 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 1;
  /* Now define structured type and commit it */
  MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &particletype);
  MPI_Type_commit(&particletype);

  compute_neighbors();   // Needed for shipping
  init_seed_times();


  // load_particle_positions(nc, last_record);  // IF WE ARE RESTARTING

  std::string input_file;
  int start = 0;
  bool boot = false;
  bool use_input_file = find_pism_input(input_file, boot, start);

  options::String poo("-o", "Output file name", "unnamed.nc");
  tracer_log_created = poo;
  tracer_log_deleted = poo;
  tracer_log_created.replace(tracer_log_created.end()-3,tracer_log_created.end(), "_tracers_created.nc");
  tracer_log_deleted.replace(tracer_log_deleted.end()-3,tracer_log_deleted.end(), "_tracers_deleted.nc");
  PIO tracer_log_nc (m_grid->com, m_config->get_string("output_format"));
  tracer_log_nc.open(tracer_log_created, PISM_READWRITE_MOVE);
  prepare_tracer_log_file(tracer_log_nc);
  tracer_log_nc.close();
  tracer_log_nc.open(tracer_log_deleted, PISM_READWRITE_MOVE);
  prepare_tracer_log_file(tracer_log_nc);
  tracer_log_nc.close();

  options::String sm("-seed_mask_file", "File for seed mask, variable name is seed_mask\n 0 = don't seed, n > 0 = seed n tracers per round.");
  m_seed_mask.create(m_grid, "seed_mask", WITHOUT_GHOSTS);
  m_seed_mask.set_attrs("internal",
                 "particles to seed in grid cell",
                 "", "");
  if (sm.is_set()){
    m_log->message(2,
		   "LAGRANGE: Reading seed_mask from file %s", sm.value().c_str());
    m_seed_mask.regrid(sm.value(), CRITICAL);
  } else {
    m_log->message(2,
		   "LAGRANGE: Setting seed_mask to one\n");
      m_seed_mask.set(1);
  }

  bool bootstrap_tracers = options::Bool("-bootstrap_tracers", "Load tracers from file while bootstraping");
  if (use_input_file && ((!boot) || bootstrap_tracers) ) {
    PIO nc(m_grid->com, "guess_mode");
    nc.open(input_file, PISM_READONLY);
    load_particle_positions(input_file);
  } else {
  tracer_counter = 0; 
  //set_initial_distribution_parameters();
  }

}
  

/** Create an initial state (cloud of "particles").
 *
 * Use this for testing or bootstrapping only.
 * Deactivated in defauld code.
 * 
 */
void PISMLagrange::set_initial_distribution_parameters() {

  verbPrintf(2, m_grid->com,
                    "  Setting initial distribution parameters...\n"); 
  unsigned int n_particles = (unsigned int)  (m_grid->Mz()-1) * (m_grid->xm()) * (m_grid->ym());
  

  // Initialize particle positions and corresponding enhancement factors
  particles.resize(n_particles);
  int i = get_offset(n_particles); // Globally coordinated numbering scheme
  if (m_grid->rank() == 0){
    tracer_counter = i; 
      i=0 ; // NO OFFSET FOR FIRST PROCESS
  }

  for (unsigned int xx = m_grid->xs(); xx < (unsigned int)m_grid->xs()+m_grid->xm(); ++xx) {
    for (unsigned int yy = m_grid->ys(); yy < (unsigned int)m_grid->ys()+m_grid->ym(); ++yy) {
      for (unsigned int zz = 0; zz < m_grid->Mz()-1; ++zz) {
	static std::list<Particle>::iterator it = particles.begin();
	  if (it==particles.end())
	  break;
	it->x = m_grid->x(xx);
	it->y = m_grid->y(yy);
	it->z = (m_grid->z(zz)+m_grid->z(zz+1))/2.0;
	it->id = i;
	it++;
	i++;
      }
    }
  }
  log_tracers(particles.begin(), particles.size(), m_grid->ctx()->time()->current(), tracer_log_created);

}

  /** Compute ranks of neighboring processes. 
   *    Necessary for passing tracers from one process to the other.
   */  
void PISMLagrange::compute_neighbors(){

  PetscInt dim, M, N, P, m, n, p, dof, s;
  DMBoundaryType bx, by, bz;
  DMDAStencilType st;
  petsc::DM::Ptr dm = m_grid->get_dm(1,0);
  DMDAGetInfo(* dm , &dim, &M, &N, &P, &m, &n, &p, &dof, &s, &bx, &by, &bz, &st);
  int r = m_grid->rank();
  int rp =m ,  rm = m  ; 
  int cp = 1 , cm = 1  ; 
  
  if (r / m +1 >=n)
    rp -= m*n;
  if (r / m ==0 )
    rm -=m*n;
  if (r%m == 0 )
      cm -= m;
  if (r%m == m-1 )
    cp -= m;


  neighbors [0] = r-rm-cm;
  neighbors [1] = r-cm;
  neighbors [2] = r+rp-cm;
  
  neighbors [3] = r-rm;
  neighbors [4] = r;
  neighbors [5] = r+rp;
  
  neighbors [6] = r-rm+cp;
  neighbors [7] = r+cp;
  neighbors [8] = r+rp+cp;

  //  for (int i=0 ; i< 9 ; i++)   std::cout <<" for "<< m_grid->rank() << "neigbor "<< i << " has rank "<< neighbors[i]<< "\n";
}

  /** Ship tracers from one process to the other 
   *
   */
  void PISMLagrange::ship_tracers(){
    std::vector<std::vector<Particle> > ship(9);
      
    for (std::list<Particle>::iterator it = particles.begin() ; it != particles.end() ;){
      const unsigned int dest =  whereto(it->x, it->y,it->z);
      if (dest != 4){
	ship[dest].push_back(*it);
	it = particles.erase(it);
      } else
	it++;
    }
    int num_ship[9];
    int num_rec[9];
    MPI_Status status[9];
    MPI_Request sreq[9], rreq[9];
    std::vector<Particle>received_tracers[9];
    for (int i = 0 ; i< 9 ; i++){
      if (i==4) continue; // not shipping to myself
      num_ship[i] = ship[i].size();
      MPI_Isend ( &num_ship[i], 1, MPI_INT, neighbors[i],   i, m_grid->com, &sreq[i]);
      MPI_Irecv ( &num_rec[i], 1, MPI_INT, neighbors[i], 8-i, m_grid->com , &rreq[i] );
    }
    for (int i = 0 ; i< 9 ; i++){
      if (i==4) continue; // not shipping to myself
      MPI_Wait(&sreq[i], &status[i]);
      MPI_Wait(&rreq[i], &status[i]);
    }

    for (int i = 0 ; i< 9 ; i++){
      if (i==4) continue; // not shipping to myself
      received_tracers[i].resize(num_rec[i]);
      MPI_Isend ( &ship[i][0], num_ship[i], particletype, neighbors[i],  i, m_grid->com, &sreq[i]);
      MPI_Irecv ( &received_tracers[i][0], num_rec[i], particletype, neighbors[i], 8-i, m_grid->com, &rreq[i] );
    }

    for (int i = 0 ; i< 9 ; i++){
      if (i==4) continue; // not shipping to myself
      MPI_Wait(&sreq[i], &status[i]);
      MPI_Wait(&rreq[i], &status[i]);
      for (std::vector<Particle>::iterator it = received_tracers[i].begin(); it != received_tracers[i].end(); it++)
	particles.push_back(*it);

      received_tracers[i].resize(0);
      ship[i].resize(0);
    }
    
  }
  
  /** To which neighbor should a particle go? 
   *   6 7 8 
   *   3 4 5 
   *   0 1 2 
   *   4 is stay local vert = y, hor = x
  */
  
  int PISMLagrange::whereto(double x, double y, double z){
    int row = 1, column = 1;
    int xs=m_grid->xs(), xm=m_grid->xm(),
      ys=m_grid->ys(), ym=m_grid->ym();
    double dx = m_grid->dx(), dy = m_grid->dy();
    if (x < m_grid->x(xs)-dx/2.)
      column=0;
    else if (x > m_grid->x(xs)+dx*(xm-0.5))
      column=2;
    
    if (y < m_grid->y(ys)-dy/2.)
      row=0;
    else if (y > m_grid->y(ys)+dy*(ym-0.5))
      row=2;

    return row*3+column;
  }
  

  void PISMLagrange::save_diagnostics(const PIO &nc) {

  }


  void PISMLagrange::save_particle_positions(const PIO &nc) {
  unsigned int n_records   = nc.inq_nrecords(); 
  unsigned int last_record = n_records - 1;
  if (not nc.inq_dimlen("tracer_id")){
    m_log->message(2,
		   "Empty tracer dimension - nothing to write");
    return;
  }
  unsigned int start, count;
  start = 0;
  count = particles.size();

  nc.put_att_double("PISM_GLOBAL","tracer_counter", PISM_INT, tracer_counter);
  // Will let rank 0 do the put, so the correct number is stored.

  
  start = get_offset(count);
  if (m_grid->rank() == 0)
    start = 0 ;
  
  std::vector <double>
    m_tracer_x(particles.size()),
    m_tracer_y(particles.size()),
    m_tracer_z(particles.size()),
    m_tracer_id(particles.size());
  std::vector<double>::iterator
    ix=m_tracer_x.begin(),
    iy=m_tracer_y.begin(),
    iz=m_tracer_z.begin(),
    id=m_tracer_id.begin();
  for (std::list<Particle>::iterator it = particles.begin() ; it != particles.end() ;it++)
    {
      * ix++ = it->x ;
      * iy++ = it->y ;
      * iz++ = it->z ;
      * id++ = (double) it->id ;
    }

    nc.put_1d_var("tracer_x", start, count, m_tracer_x); 
    nc.put_1d_var("tracer_y", start, count, m_tracer_y); 
    nc.put_1d_var("tracer_z", start, count, m_tracer_z);
    nc.put_1d_var("tracer_id", start, count, m_tracer_id); 
  }
  
  void PISMLagrange::load_particle_positions(const std::string input_file) {
    PIO nc(m_grid->com, "guess_mode");
    nc.open(input_file, PISM_READONLY);


    unsigned int n_records = nc.inq_nrecords();
    unsigned int last_record = n_records - 1;

    std::vector<double> num_trac(1, 0.0);
    
    tracer_counter = (unsigned int) nc.get_att_double("PISM_GLOBAL","tracer_counter")[0];
    unsigned int start, count;
    start = 0;
    count = nc.inq_dimlen("tracer_id");

    // OPTION FOR MANY TRACERS: READ up to  10/100/... Million at a time and sort them
    // before reading the next batch. -- ::TODO::
    
    const unsigned int n = count;
    std::vector<double> x(n), y(n),z(n), id(n);

    
    nc.get_1d_var("tracer_x", start, count, x);
    nc.get_1d_var("tracer_y", start, count, y);
    nc.get_1d_var("tracer_z", start, count, z); 
    nc.get_1d_var("tracer_id", start, count, id); 
    
    std::vector<double>::iterator
      ix = x.begin(),
      iy = y.begin(),
      iz = z.begin(),
      iid = id.begin();
    while (ix != x.end()){
      if (whereto(*ix, *iy, *iz) == 4 ){
	if (*iz < 0 ) continue;
	Particle p = {*ix, *iy, *iz, (int) *iid};
	particles.push_back(p);
      }
      ix++; iy++ ; iz++ ; iid++;
    }
    
  }
  

  /** coordinate offsets in tracer numbering - FUNKY! READ DESCRIPTION
   *
   * @return RANK=0 receives the total, all others get the sum of the prev.
   *
   *
   */
  
  int PISMLagrange::get_offset(const int contribution) {

    int retval = 0;

    if (m_grid->rank()==0){ // Main thread has to handle the mess.

      //Init MPI stuff
      int mpi_comm_size = 0 ;
      MPI_Comm_size(m_grid->com, &mpi_comm_size);
      MPI_Request req[mpi_comm_size];
      MPI_Status status[mpi_comm_size];

      // Fields for counting and such
      int counts[mpi_comm_size], offsets[mpi_comm_size];
      offsets[0] = 0;
       counts[0] = contribution;

      for (int i = 1 ; i < mpi_comm_size; i++){
	MPI_Irecv ( &counts[i], 1, MPI_INT, i, 0, m_grid->com, &req[i]);
      } // Just collecting numbers. Can happen in parallel.

      for (int i = 1 ; i < mpi_comm_size; i++){
	MPI_Wait(&req[i], &status[i]);
	offsets[i] = offsets[i-1] + counts[i-1];
	MPI_Isend ( &offsets[i], 1, MPI_INT, i, 0, m_grid->com, &req[i]);
      } // We only need this and all previos receives completed.
      

      for (int i = 1 ; i < mpi_comm_size; i++){
	MPI_Wait ( &req[i], &status[i]);
      } // Wait till everything is shipped
      retval = offsets[mpi_comm_size-1]+counts[mpi_comm_size-1];
      return retval;  // ON RANK ZERO WE RETURN THE TOTAL!
      
    } else {
      MPI_Status status;
      MPI_Request req;
      MPI_Isend ( &contribution, 1, MPI_INT, 0, 0, m_grid->com, &req);
      MPI_Wait ( &req, &status);
      MPI_Irecv ( &retval, 1, MPI_INT, 0, 0, m_grid->com, &req);
      MPI_Wait ( &req, &status);
    }
  
  return retval;
  }

  /** processes -seed_times argument analogue to -extra_times. Extra: -seed_seasons
   * 
   * -seed_seasons will generate 4 extra seeds for 3, 6, 9, 12 months after original seeding.
   */ 
  void PISMLagrange::init_seed_times(){
    last_seed = -9e20;               // will be set in write_extras()
    options::String times("-seed_times", "Specifies times to save at");
    bool seed_seasons =  options::Bool("-seed_seasons", "Seed 4 seasons (.0, .25, .5, .75 for each seeding time");
    

    if (times.is_set()){
      try {
	m_grid->ctx()->time()->parse_times(times, seed_times);
      } catch (RuntimeError &e) {
	e.add_context("parsing the -seed_times argument %s", times->c_str());
	throw;
      }
      
      if (seed_times.size() == 0) {
	throw RuntimeError("-seed_times set, but no seeding during the run specified.");
      }
      
      
    }
    if (seed_seasons){
      std::vector<double> new_seed_times(5* seed_times.size());
      std::vector<double>::iterator nst = new_seed_times.begin();
      for (std::vector<double>::iterator it = seed_times.begin() ; it != seed_times.end() ; it++)
	{
	  const double y = m_grid->ctx()->time()->increment_date(*it, 1)- *it;
	  *(nst++) = *it;
	  *(nst++) = *it + .25 * y;
	  *(nst++) = *it + .5  * y;
	  *(nst++) = *it + .75 * y;
	  *(nst++) = *it + y;
	}
      seed_times=new_seed_times;
    }
    next_seed = seed_times.begin();
    while (next_seed != seed_times.end() && *next_seed < m_grid->ctx()->time()->start()-1)
      {
	next_seed++;}

  }


  bool PISMLagrange::check_seed(){
    double seed_after = -1.0e30;
    std::vector<double>::iterator current_seed;
    if (next_seed != seed_times.end() &&
	(*next_seed <= m_grid->ctx()->time()->current() ||
	 fabs(m_grid->ctx()->time()->current() - *next_seed) < 1.0)) {
    // the condition above is "true" if we passed a requested time or got to
    // within 1 second from it
      
      current_seed = next_seed;
      
      // update next_extra
      while (next_seed != seed_times.end() &&
	     (*next_seed <= m_grid->ctx()->time()->current() ||
            fabs(m_grid->ctx()->time()->current() - *next_seed) < 1.0)) {
	next_seed++;
      }
      
      seed_after = *current_seed;
    } else {
      return false;
    }
    
    if (seed_after < m_grid->ctx()->time()->start()) {
      // Suppose a user tells PISM to write data at times 0:1000:10000. Suppose
      // also that PISM writes a backup file at year 2500 and gets stopped.
      //
      // When restarted, PISM will decide that it's time to write data for time
      // 2000, but
      // * that record was written already and
      // * PISM will end up writing at year 2500, producing a file containing one
      //   more record than necessary.
      //
      // This check makes sure that this never happens.
      return false;
    }
    return true; 
  }

  /** Seed tracers on ice surface
   *
   * Puts tracers on ice surface at current time t (not t+dt). 
   * One tracer per grid point centered in the cell.
   */
  void PISMLagrange::seed(){

    const Profiling &profiling = m_grid->ctx()->profiling();
    profiling.begin("tracer seeding");
    std::mt19937 rng;

    const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
    IceModelVec::AccessList list;
    list.add(*thickness);
    list.add(m_seed_mask);
    
    std::list<Particle> new_tracers;
    int id = 0 ;
    for (Points p(*m_grid); p; p.next()){
      int count = 0 ;
      const int i = p.i(), j = p.j();
      const double thk = (*thickness)(i,j);
      while (thk > 0 && m_seed_mask(i,j) > count ){
	double ox = 0, oy = 0 ;
	if (count > 0 ){
	  if (count == 1){
	    rng.seed(i*m_grid->My() + j);
	  }
	  ox = ((rng() / (double) (rng.max())) -.5) * m_grid->dx();
	  oy = ((rng() / (double) (rng.max())) -.5) * m_grid->dy();
	}
	const double xx = m_grid->x(i) + ox,
                     yy = m_grid->y(j) + oy;
	Particle part = {xx , yy , evaluate_at_point(*thickness, xx, yy), id };
	new_tracers.push_back(part);
	id++;
	count++;
	}
    }
    int offset;
    if (m_grid->rank() == 0){
      offset = tracer_counter ;
      tracer_counter = get_offset( offset + id);
    } else{
      offset = get_offset(id);
    }
    for (std::list<Particle>::iterator it = new_tracers.begin() ; it != new_tracers.end() ;it++)
      it->id+=offset;
    particles.insert(particles.end(), new_tracers.begin(), new_tracers.end());

    log_tracers(new_tracers.begin(), new_tracers.size(), m_grid->ctx()->time()->current(), tracer_log_created);
    
    profiling.end("tracer seeding");
  
  }

  /** Any tracer that's more than 1 mm above its cell's ice surface is discarded.
   * no interpolation. Pure cell center.
   * 
   */
  
  void PISMLagrange::remove_flying_tracers(){

    const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
    IceModelVec::AccessList list;
    list.add(*thickness);

    std::list<Particle> deleted;

    const double
      x0 = m_grid->x(0),
      y0 = m_grid->y(0),
      dx = m_grid->dx(),
      dy = m_grid->dy();
    for (std::list<Particle>::iterator it = particles.begin() ; it != particles.end() ;)
      {
	if (it->z > evaluate_at_point(*thickness, it->x,it->y) +1e-3){ // 1 mm grace, let's call it surfce roughness ;)
	  deleted.push_back(*it);
	  it = particles.erase(it);
	} else
	  it++;
      }
    log_tracers(deleted.begin(), deleted.size(), m_grid->ctx()->time()->current(), tracer_log_deleted);
  }

  
  void PISMLagrange::remove_bedrock_tracers(){


    std::list<Particle> deleted;
    
    for (std::list<Particle>::iterator it = particles.begin() ; it != particles.end() ;it++)
      {      
	if (it->z < 1e-6){ // 1 micrometer above ground is in ground
	  deleted.push_back(*it);
	  it = particles.erase(it);
	  it-- ;
	}
      }
    log_tracers(deleted.begin(), deleted.size(), m_grid->ctx()->time()->current(), tracer_log_deleted);
  }

  void PISMLagrange::prepare_tracer_log_file(const PIO & nc){

    std::vector<double> index;
    bool dim_exists = false;
    dim_exists = nc.inq_dim("id");
    if (dim_exists == false) {
      nc.redef();
      nc.def_dim("tracer_id", PISM_UNLIMITED );
    }


    std::vector<std::string> dims(1);
    dims[0] = "tracer_id";
    // particle positions
    std::vector<std::string> fields ;
    fields.push_back("tracer_x");
    fields.push_back("tracer_y");
    fields.push_back("tracer_z");
    fields.push_back("tracer_id");
    fields.push_back("time");

    for (std::vector<std::string>::iterator it = fields.begin(); it != fields.end(); ++it) {
      if (not nc.inq_var(*it)) {
	nc.redef();
	nc.def_var(*it, PISM_DOUBLE, dims);
      }
    }


  }

void PISMLagrange::log_tracers(const std::list<Particle>::iterator first , const size_t count, const double a_time, const std::string filename ){

    unsigned int my_offset = get_offset(count); // Will be incremented by File offset later.
    if (m_grid->rank() == 0 )
      my_offset = 0;

    std::vector<double>
      x(count),
      y(count),
      z(count),
      id(count),
      time(count, a_time);
    std::list<Particle>::iterator  it = first;
    for (size_t i = 0 ; i < count ; i++)
      {
	x[i]=it->x; 
	y[i]=it->y; 
	z[i]=it->z; 
	id[i]=(double) it->id;
	it++;
      }
    
    
    PIO nc (m_grid->com, "guess_mode");
    nc.open(filename, PISM_READWRITE);
    my_offset += nc.inq_dimlen("tracer_id");
    nc.put_1d_var("tracer_id", my_offset, count, id);
    nc.put_1d_var("tracer_x", my_offset, count, x);
    nc.put_1d_var("tracer_y", my_offset, count, y);
    nc.put_1d_var("tracer_z", my_offset, count, z);
    nc.put_1d_var("time", my_offset, count, time);
  }

  void PISMLagrange::print_tracer(const Particle & p){
    std::cout <<" Particle " << p.id << "\t x="<< p.x << "\t y="<< p.y << "\t z="<< p.z<<"\n";
  }
} // end of namespace pism
