// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev and Ed Bueler
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

#include "ShallowStressBalance.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/rheology/FlowLawFactory.hh"

#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "base/util/IceModelVec2CellType.hh"

#include "SSB_diagnostics.hh"

namespace pism {
namespace stressbalance {

using pism::mask::ice_free;

ShallowStressBalance::ShallowStressBalance(IceGrid::ConstPtr g)
  : Component(g), m_basal_sliding_law(NULL), m_flow_law(NULL), m_EC(g->ctx()->enthalpy_converter()) {

  m_bc_values = NULL;
  m_bc_mask = NULL;
  m_sea_level = 0.0;

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  if (m_config->get_boolean("basal_resistance.pseudo_plastic.enabled") == true) {
    m_basal_sliding_law = new IceBasalResistancePseudoPlasticLaw(*m_config);
  } else {
    m_basal_sliding_law = new IceBasalResistancePlasticLaw(*m_config);
  }

  m_velocity.create(m_grid, "bar", WITH_GHOSTS, WIDE_STENCIL); // components ubar, vbar
  m_velocity.set_attrs("model_state",
                       "thickness-advective ice velocity (x-component)", 
                       "m s-1", "", 0);
  m_velocity.set_attrs("model_state",
                       "thickness-advective ice velocity (y-component)",
                       "m s-1", "", 1);

  m_velocity.metadata(0).set_string("glaciological_units", "m year-1");
  m_velocity.metadata(1).set_string("glaciological_units", "m year-1");

  m_velocity.write_in_glaciological_units = true;

  m_basal_frictional_heating.create(m_grid, "bfrict", WITHOUT_GHOSTS);
  m_basal_frictional_heating.set_attrs("diagnostic",
                                       "basal frictional heating",
                                       "W m-2", "");
  m_basal_frictional_heating.metadata().set_string("glaciological_units", "mW m-2");
  m_basal_frictional_heating.write_in_glaciological_units = true;
}

ShallowStressBalance::~ShallowStressBalance() {
  delete m_basal_sliding_law;
}

void ShallowStressBalance::init() {
  this->init_impl();
}

void ShallowStressBalance::init_impl() {
  // empty
}

std::string ShallowStressBalance::stdout_report() const {
  return "";
}

const rheology::FlowLaw* ShallowStressBalance::flow_law() const {
  return m_flow_law;
}

EnthalpyConverter::Ptr ShallowStressBalance::enthalpy_converter() const {
  return m_EC;
}

const IceBasalResistancePlasticLaw* ShallowStressBalance::sliding_law() const {
  return m_basal_sliding_law;
}


void ShallowStressBalance::set_boundary_conditions(const IceModelVec2Int &locations,
                                                   const IceModelVec2V &velocities) {
  m_bc_values = &velocities;
  m_bc_mask = &locations;
}

//! \brief Get the thickness-advective 2D velocity.
const IceModelVec2V& ShallowStressBalance::velocity() const {
  return m_velocity;
}

//! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
const IceModelVec2S& ShallowStressBalance::basal_frictional_heating() {
  return m_basal_frictional_heating;
}


std::map<std::string, Diagnostic::Ptr> ShallowStressBalance::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"beta",     Diagnostic::Ptr(new SSB_beta(this))},
    {"taub",     Diagnostic::Ptr(new SSB_taub(this))},
    {"taub_mag", Diagnostic::Ptr(new SSB_taub_mag(this))},
    {"taud",     Diagnostic::Ptr(new SSB_taud(this))},
    {"taud_mag", Diagnostic::Ptr(new SSB_taud_mag(this))}
  };
  return result;
}


ZeroSliding::ZeroSliding(IceGrid::ConstPtr g)
  : ShallowStressBalance(g) {

  // Use the SIA flow law.
  rheology::FlowLawFactory ice_factory("stress_balance.sia.", m_config, m_EC);
  m_flow_law = ice_factory.create();
}

ZeroSliding::~ZeroSliding() {
  delete m_flow_law;
}

//! \brief Update the trivial shallow stress balance object.
void ZeroSliding::update(bool fast, double sea_level, const IceModelVec2S &melange_back_pressure) {
  (void) sea_level;
  (void) melange_back_pressure;

  if (not fast) {
    m_velocity.set(0.0);
    m_basal_frictional_heating.set(0.0);
  }
}

//! \brief Compute the basal frictional heating.
/*!
  Ice shelves have zero basal friction heating.

  \param[in] V *basal* sliding velocity
  \param[in] tauc basal yield stress
  \param[in] mask (used to determine if floating or grounded)
  \param[out] result
 */
void ShallowStressBalance::compute_basal_frictional_heating(const IceModelVec2V &V,
                                                            const IceModelVec2S &tauc,
                                                            const IceModelVec2CellType &mask,
                                                            IceModelVec2S &result) const {

  IceModelVec::AccessList list{&V, &result, &tauc, &mask};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i,j)) {
      result(i,j) = 0.0;
    } else {
      const double
        C = m_basal_sliding_law->drag(tauc(i,j), V(i,j).u, V(i,j).v),
        basal_stress_x = - C * V(i,j).u,
        basal_stress_y = - C * V(i,j).v;
      result(i,j) = - basal_stress_x * V(i,j).u - basal_stress_y * V(i,j).v;
    }
  }
}


//! \brief Compute 2D deviatoric stresses.
/*! Note: IceModelVec2 result has to have dof == 3. */
void ShallowStressBalance::compute_2D_stresses(const IceModelVec2V &V,
                                               const IceModelVec2CellType &mask,
                                               IceModelVec2 &result) const {
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  if (result.get_ndof() != 3) {
    throw RuntimeError(PISM_ERROR_LOCATION, "result.get_dof() == 3 is required");
  }

  // NB: uses constant ice hardness; choice is to use SSA's exponent; see issue #285
  double hardness = pow(m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
                        -1.0 / m_config->get_double("stress_balance.ssa.Glen_exponent"));

  IceModelVec::AccessList list{&V, &result, &mask};

  //////////////////////////////////////////////////////////////////////////////////////
  bool mean_hardness = options::Bool("-use_mean_hardness", "Use mean hardness for deviatoric stress calcuation");

  //if (mean_hardness){

    const IceModelVec3  *enthalpy  = m_grid->variables().get_3d_scalar("enthalpy");
    const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

    list.add(*enthalpy);
    list.add(*thickness);

    const double *z = &m_grid->z()[0];
  //}
  /////////////////////////////////////////////////////////////////////////////////////


  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_free(mask.as_int(i,j))) {
      result(i,j,0) = 0.0;
      result(i,j,1) = 0.0;
      result(i,j,2) = 0.0;
      continue;
    }

    StarStencil<int> m = mask.int_star(i,j);
    StarStencil<Vector2> U = V.star(i,j);

    // strain in units s-1
    double u_x = 0, u_y = 0, v_x = 0, v_y = 0,
      east = 1, west = 1, south = 1, north = 1;

    // Computes u_x using second-order centered finite differences written as
    // weighted sums of first-order one-sided finite differences.
    //
    // Given the cell layout
    // *----n----*
    // |         |
    // |         |
    // w         e
    // |         |
    // |         |
    // *----s----*
    // east == 0 if the east neighbor of the current cell is ice-free. In
    // this case we use the left- (west-) sided difference.
    //
    // If both neighbors in the east-west (x) direction are ice-free the
    // x-derivative is set to zero (see u_x, v_x initialization above).
    //
    // Similarly in y-direction.
    if (ice_free(m.e)) {
      east = 0;
    }
    if (ice_free(m.w)) {
      west = 0;
    }
    if (ice_free(m.n)) {
      north = 0;
    }
    if (ice_free(m.s)) {
      south = 0;
    }

    if (west + east > 0) {
      u_x = 1.0 / (dx * (west + east)) * (west * (U.ij.u - U[West].u) + east * (U[East].u - U.ij.u));
      v_x = 1.0 / (dx * (west + east)) * (west * (U.ij.v - U[West].v) + east * (U[East].v - U.ij.v));
    }

    if (south + north > 0) {
      u_y = 1.0 / (dy * (south + north)) * (south * (U.ij.u - U[South].u) + north * (U[North].u - U.ij.u));
      v_y = 1.0 / (dy * (south + north)) * (south * (U.ij.v - U[South].v) + north * (U[North].v - U.ij.v));
    }

    //////////////////////////////////////////////////////////////////
    if (mean_hardness) {
      double H = (*thickness)(i, j);
      unsigned int k = m_grid->kBelowHeight(H);
      hardness = averaged_hardness(*m_flow_law, H, k, &z[0], enthalpy->get_column(i, j));
    }
    ////////////////////////////////////////////////////////////////////

    double nu = 0.0;
    m_flow_law->effective_viscosity(hardness,
                                    secondInvariant_2D(Vector2(u_x, v_x), Vector2(u_y, v_y)),
                                    &nu, NULL);

    //get deviatoric stresses
    result(i,j,0) = nu*u_x;
    result(i,j,1) = nu*v_y;
    result(i,j,2) = 0.5*nu*(u_y+v_x);
  }
}

SSB_taud::SSB_taud(const ShallowStressBalance *m)
  : Diag<ShallowStressBalance>(m) {

  m_dof = 2;

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "taud_x"),
            SpatialVariableMetadata(m_sys, "taud_y")};

  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < m_dof; ++k) {
    m_vars[k].set_string("comment",
                       "this field is purely diagnostic (not used by the model)");
  }
}

/*!
 * The driving stress computed here is not used by the model, so this
 * implementation intentionally does not use the eta-transformation or special
 * cases at ice margins.
 */
IceModelVec::Ptr SSB_taud::compute_impl() {

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "result", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2S *surface = m_grid->variables().get_2d_scalar("surface_altitude");

  double standard_gravity = m_config->get_double("constants.standard_gravity"),
    ice_density = m_config->get_double("constants.ice.density");

  IceModelVec::AccessList list{surface, thickness, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double pressure = ice_density * standard_gravity * (*thickness)(i,j);
    if (pressure <= 0.0) {
      (*result)(i,j).u = 0.0;
      (*result)(i,j).v = 0.0;
    } else {
      (*result)(i,j).u = - pressure * surface->diff_x_p(i,j);
      (*result)(i,j).v = - pressure * surface->diff_y_p(i,j);
    }
  }

  return result;
}

SSB_taud_mag::SSB_taud_mag(const ShallowStressBalance *m)
  : Diag<ShallowStressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "taud_mag")};

  set_attrs("magnitude of the gravitational driving stress at the base of ice", "",
            "Pa", "Pa", 0);
  m_vars[0].set_string("comment",
                     "this field is purely diagnostic (not used by the model)");
}

IceModelVec::Ptr SSB_taud_mag::compute_impl() {

  // Allocate memory:
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "taud_mag", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec2V::Ptr taud = IceModelVec2V::ToVector(SSB_taud(model).compute());

  result->set_to_magnitude(*taud);

  return result;
}

SSB_taub::SSB_taub(const ShallowStressBalance *m)
  : Diag<ShallowStressBalance>(m) {
  m_dof = 2;

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "taub_x"),
            SpatialVariableMetadata(m_sys, "taub_y")};

  set_attrs("X-component of the shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < m_dof; ++k) {
    m_vars[k].set_string("comment",
                       "this field is purely diagnostic (not used by the model)");
  }
}


IceModelVec::Ptr SSB_taub::compute_impl() {

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "result", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->metadata(1) = m_vars[1];

  const IceModelVec2V        &velocity = model->velocity();
  const IceModelVec2S        *tauc     = m_grid->variables().get_2d_scalar("tauc");
  const IceModelVec2CellType &mask     = *m_grid->variables().get_2d_cell_type("mask");

  const IceBasalResistancePlasticLaw *basal_sliding_law = model->sliding_law();

  IceModelVec::AccessList list{tauc, &velocity, &mask, result.get()};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.grounded_ice(i,j)) {
      double beta = basal_sliding_law->drag((*tauc)(i,j), velocity(i,j).u, velocity(i,j).v);
      (*result)(i,j).u = - beta * velocity(i,j).u;
      (*result)(i,j).v = - beta * velocity(i,j).v;
    } else {
      (*result)(i,j).u = 0.0;
      (*result)(i,j).v = 0.0;
    }
  }

  return result;
}

SSB_taub_mag::SSB_taub_mag(const ShallowStressBalance *m)
  : Diag<ShallowStressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "taub_mag")};

  set_attrs("magnitude of the basal shear stress at the base of ice",
            "magnitude_of_land_ice_basal_drag", // InitMIP "standard" name
            "Pa", "Pa", 0);
  m_vars[0].set_string("comment",
                     "this field is purely diagnostic (not used by the model)");
}

IceModelVec::Ptr SSB_taub_mag::compute_impl() {

  // Allocate memory:
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "taub_mag", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec2V::Ptr taub = IceModelVec2V::ToVector(SSB_taub(model).compute());

  result->set_to_magnitude(*taub);

  return result;
}

/**
 * Shallow stress balance class that reads `u` and `v` fields from a
 * file and holds them constant.
 *
 * The only use I can think of right now is testing.
 */
PrescribedSliding::PrescribedSliding(IceGrid::ConstPtr g)
  : ZeroSliding(g) {
  // empty
}

PrescribedSliding::~PrescribedSliding() {
  // empty
}

void PrescribedSliding::update(bool fast, double sea_level,
                               const IceModelVec2S &melange_back_pressure) {
  (void) sea_level;
  (void) melange_back_pressure;
  if (not fast) {
    m_basal_frictional_heating.set(0.0);
  }
}

void PrescribedSliding::init_impl() {
  ShallowStressBalance::init_impl();

  options::String input_filename("-prescribed_sliding_file",
                                 "name of the file to read velocity fields from");
  if (not input_filename.is_set()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "option -prescribed_sliding_file is required.");
  }

  m_velocity.regrid(input_filename, CRITICAL);
}

SSB_beta::SSB_beta(const ShallowStressBalance *m)
  : Diag<ShallowStressBalance>(m) {
  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "beta")};

  set_attrs("basal drag coefficient", "", "Pa s / m", "Pa s / m", 0);
}

IceModelVec::Ptr SSB_beta::compute_impl() {

  // Allocate memory:
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "beta", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->write_in_glaciological_units = true;

  const IceModelVec2S *tauc = m_grid->variables().get_2d_scalar("tauc");

  const IceBasalResistancePlasticLaw *basal_sliding_law = model->sliding_law();

  const IceModelVec2V &velocity = model->velocity();

  IceModelVec::AccessList list{tauc, &velocity, result.get()};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) =  basal_sliding_law->drag((*tauc)(i,j), velocity(i,j).u, velocity(i,j).v);
  }

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
