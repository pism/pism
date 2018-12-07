#include <iostream>

#include <pism/energy/BedThermalUnit.hh>
#include <pism/util/EnthalpyConverter.hh>
#include <pism/util/io/PIO.hh>
#include <pism/util/io/io_helpers.hh>
#include "pism/energy/EnergyModel.hh"

#include "pism/icebin/IBIceModel.hh"
#include "pism/icebin/IBSurfaceModel.hh"

#include "pism/coupler/SeaLevel.hh"
#include "pism/coupler/ocean/sea_level/Initialization.hh"

namespace pism {
namespace icebin {

// ================================

IBIceModel::IBIceModel(IceGrid::Ptr g, Context::Ptr context, IBIceModel::Params const &_params)
    : pism::IceModel(g, context), params(_params) {
  // empty
}

IBIceModel::~IBIceModel() {
  // empty
}


void IBIceModel::allocate_subglacial_hydrology() {
  printf("BEGIN IBIceModel::allocate_subglacial_hydrology()\n");

  if (pism::IceModel::m_subglacial_hydrology) {
    return; // indicates it has already been allocated
  }

  m_subglacial_hydrology.reset(new pism::hydrology::NullTransport(m_grid));

  m_submodels["subglacial hydrology"] = m_subglacial_hydrology.get();

  printf("END IBIceModel::allocate_subglacial_hydrology()\n");
}


void IBIceModel::allocate_couplers() {
  // Initialize boundary models:

  if (not m_surface) {

    m_log->message(2, "# Allocating a surface process model or coupler...\n");

    m_surface.reset(new IBSurfaceModel(m_grid));

    m_submodels["surface process model"] = m_surface.get();
  }

  if (not m_ocean) {
    m_log->message(2, "# Allocating an ocean model or coupler...\n");

    ocean::Factory po(m_grid);
    m_ocean = po.create();

    m_submodels["ocean model"] = m_ocean.get();
  }

  if (not m_sea_level) {
    using namespace ocean::sea_level;
    std::shared_ptr<SeaLevel> sea_level(new SeaLevel(m_grid));
    m_sea_level.reset(new InitializationHelper(m_grid, sea_level));
    m_submodels["sea level forcing"] = m_sea_level.get();
  }
}


void IBIceModel::allocate_storage() {
  super::allocate_storage();

  printf("BEGIN IBIceModel::allocate_storage()\n");
  base.create(m_grid, "", WITHOUT_GHOSTS);
  cur.create(m_grid, "", WITHOUT_GHOSTS);
  rate.create(m_grid, "", WITHOUT_GHOSTS);
  printf("END IBIceModel::allocate_storage()\n");

  M1.create(m_grid, "M1", pism::WITHOUT_GHOSTS);
  M2.create(m_grid, "M2", pism::WITHOUT_GHOSTS);
  M1.create(m_grid, "H1", pism::WITHOUT_GHOSTS);
  M2.create(m_grid, "H2", pism::WITHOUT_GHOSTS);
  M1.create(m_grid, "V1", pism::WITHOUT_GHOSTS);
  M2.create(m_grid, "V2", pism::WITHOUT_GHOSTS);

  std::cout << "IBIceModel Conservation Formulas:" << std::endl;
  cur.print_formulas(std::cout);
}

void IBIceModel::energy_step() {

  printf("BEGIN IBIceModel::energy_step(t=%f, dt=%f)\n", t_TempAge, dt_TempAge);

  // Enthalpy and mass continuity are stepped with different timesteps.
  // Fish out the timestep relevant to US.
  // const double my_t0 = t_TempAge;          // Time at beginning of timestep
  const double my_dt = dt_TempAge;

  // =========== BEFORE Energy Step

  // =========== The Energy Step Itself
  super::energy_step();

  // =========== AFTER Energy Step

  // We need to integrate over strain_heating and geothermal_flux, which
  // are given in PISM as rates.


  // --------- Upward Geothermal Flux
  // Use actual geothermal flux, not the long-term average..
  // See: file:///Users/rpfische/git/pism/build/doc/browser/html/classPISMBedThermalUnit.html#details
  {
    cur.upward_geothermal_flux.add(my_dt, m_btu->flux_through_top_surface());
  }

  // ----------- Geothermal Flux
  cur.geothermal_flux.add(my_dt, m_btu->flux_through_bottom_surface());

  // ---------- Basal Frictional Heating (see iMenthalpy.cc l. 220)
  IceModelVec2S const &Rb(m_stress_balance->basal_frictional_heating());
  cur.basal_frictional_heating.add(my_dt, Rb);

  // NOTE: strain_heating is inf at the coastlines.
  // See: https://github.com/pism/pism/issues/292
  // ------------ Volumetric Strain Heating
  // strain_heating_sum += my_dt * sum_columns(strainheating3p)
  const IceModelVec3 &strain_heating3(m_stress_balance->volumetric_strain_heating());
  // cur.strain_heating = cur.strain_heating * 1.0 + my_dt * sum_columns(strain_heating3p)
  strain_heating3.sumColumns(cur.strain_heating, 1.0, my_dt);

  printf("END IBIceModel::energy_step(time=%f)\n", t_TempAge);
}

void IBIceModel::massContExplicitStep(double dt,
                                      const IceModelVec2Stag &diffusive_flux,
                                      const IceModelVec2V &advective_velocity) {

  printf("BEGIN IBIceModel::MassContExplicitStep()\n");

  _ice_density              = m_config->get_double("constants.ice.density");
  _meter_per_s_to_kg_per_m2 = dt * _ice_density;


  // =========== The Mass Continuity Step Itself
  // This will call through to accumulateFluxes_massContExplicitStep()
  // in the inner loop
  {
    AccessList access{ &cur.pism_smb,           &cur.melt_grounded, &cur.melt_floating,
                       &cur.internal_advection, &cur.href_to_h,     &cur.nonneg_rule };

    // FIXME: this is obviously broken now that PISM uses GeometryEvolution instead.
    (void) diffusive_flux;
    (void) advective_velocity;
    // super::massContExplicitStep(dt, diffusive_flux, advective_velocity);
  }

  // =========== AFTER the Mass Continuity Step

  // ----------- SMB: Pass inputs through to outputs.
  // They are needed to participate in mass/energy budget
  IBSurfaceModel *ib_surface = ib_surface_model();


  {
    AccessList access{ &ib_surface->icebin_massxfer, &ib_surface->icebin_enthxfer, &ib_surface->icebin_deltah,
                       &cur.icebin_xfer, &cur.icebin_deltah };

    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        cur.icebin_xfer.mass(i, j) += m_dt * ib_surface->icebin_massxfer(i, j);
        cur.icebin_xfer.enth(i, j) += m_dt * ib_surface->icebin_enthxfer(i, j);
        cur.icebin_deltah(i, j) += m_dt * ib_surface->icebin_deltah(i, j);
      }
    }
  }

  printf("END IBIceModel::MassContExplicitStep()\n");
}


/** This is called IMMEDIATELY after ice is gained/lost in
iMgeometry.cc (massContExplicitStep()).  Here we can record the same
values that PISM saw when moving ice around. */
void IBIceModel::accumulateFluxes_massContExplicitStep(int i, int j,
                                                       double surface_mass_balance, // [m s-1] ice equivalent
                                                       double meltrate_grounded,    // [m s-1] ice equivalent
                                                       double meltrate_floating,    // [m s-1] ice equivalent
                                                       double divQ_SIA,             // [m s-1] ice equivalent
                                                       double divQ_SSA,             // [m s-1] ice equivalent
                                                       double Href_to_H_flux,       // [m] ice equivalent
                                                       double nonneg_rule_flux)     // [m] ice equivalent
{
  EnthalpyConverter::Ptr EC = ctx()->enthalpy_converter();

  // -------------- Melting
  double p_basal             = EC->pressure(m_geometry.ice_thickness(i, j));
  double T                   = EC->melting_temperature(p_basal);
  double specific_enth_basal = EC->enthalpy_permissive(T, 1.0, p_basal);
  double mass;

  // ------- Melting at base of ice sheet
  mass = -meltrate_grounded * _meter_per_s_to_kg_per_m2;
  cur.melt_grounded.mass(i, j) += mass;
  cur.melt_grounded.enth(i, j) += mass * specific_enth_basal;

  // ------- Melting under ice shelf
  mass = -meltrate_floating * _meter_per_s_to_kg_per_m2;
  cur.melt_floating.mass(i, j) += mass;
  cur.melt_floating.enth(i, j) += mass * specific_enth_basal;


  // -------------- internal_advection
  const int ks             = m_grid->kBelowHeight(m_geometry.ice_thickness(i, j));
  // Approximate, we will use the enthalpy of the top layer...
  double specific_enth_top = m_energy_model->enthalpy().get_column(i, j)[ks];

  mass = -(divQ_SIA + divQ_SSA) * _meter_per_s_to_kg_per_m2;

  cur.internal_advection.mass(i, j) += mass;
  cur.internal_advection.enth(i, j) += mass * specific_enth_top;


  // -------------- Get the easy veriables out of the way...
  mass = surface_mass_balance * _meter_per_s_to_kg_per_m2;
  cur.pism_smb.mass(i, j) += mass;
  cur.pism_smb.enth(i, j) += mass * specific_enth_top;
  cur.nonneg_rule(i, j) -= nonneg_rule_flux * _ice_density;
  cur.href_to_h(i, j) += Href_to_H_flux * _ice_density;


  //  printf("END IBIceModel::accumulateFluxes_MassContExplicitStep()\n");
}


void IBIceModel::prepare_nc(std::string const &fname, std::unique_ptr<PIO> &nc) {

  //    nc.reset(new PIO(m_grid->com, m_grid->ctx()->config()->get_string("output.format")));

  nc.reset(new PIO(m_grid->com, m_config->get_string("output.format"),
                   fname, PISM_READWRITE_MOVE));

  io::define_time(*nc, m_grid->ctx()->config()->get_string("time.dimension_name"), m_grid->ctx()->time()->calendar(),
                  m_grid->ctx()->time()->CF_units_string(), m_grid->ctx()->unit_system());

  // These are in iMtimseries, but not listed as required in iceModelVec.hh
  //    nc->put_att_text(m_config.get_string("time.dimension_name"),
  //                           "bounds", "time_bounds");
  //    write_metadata(nc, true, false);
  //  nc->close():
}

/** @param t0 Time of last time we coupled. */
void IBIceModel::set_rate(double dt) {

  printf("BEGIN IBIceModel::set_rate(dt=%f)\n", dt);

  double by_dt = 1.0 / dt;

  compute_enth2(cur.total.enth, cur.total.mass);
  cur.set_epsilon(m_grid);

  // Compute differences, and set base = cur
  auto base_ii(base.all_vecs.begin());
  auto cur_ii(cur.all_vecs.begin());
  auto rate_ii(rate.all_vecs.begin());
  for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii, ++rate_ii) {
    IceModelVec2S &vbase(base_ii->vec);
    IceModelVec2S &vcur(cur_ii->vec);
    IceModelVec2S &vrate(rate_ii->vec);

    {
      AccessList access{ &vbase, &vcur, &vrate };
      for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
        for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
          // rate = cur - base: Just for DELTA and EPISLON flagged vectors
          if (base_ii->flags & (MassEnergyBudget::DELTA | MassEnergyBudget::EPSILON)) {
            vrate(i, j) = (vcur(i, j) - vbase(i, j)) * by_dt;
          } else {
            // Or else just copy the to rate
            vrate(i, j) = vcur(i, j);
          }
        }
      }
    }
  }

  printf("END IBIceModel::set_rate()\n");
}

void IBIceModel::reset_rate() {
  // Compute differences, and set base = cur
  auto base_ii(base.all_vecs.begin());
  auto cur_ii(cur.all_vecs.begin());
  for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii) {
    IceModelVec2S &vbase(base_ii->vec);
    IceModelVec2S &vcur(cur_ii->vec);

    // This cannot go in the loop above with PETSc because
    // vbase is needed on the RHS of the equations above.
    AccessList access{ &vbase, &vcur };
    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        // base = cur: For ALL vectors
        vbase(i, j) = vcur(i, j);
      }
    }
  }


#if 0
rate.geothermal_flux.begin_access();
printf("GG rate.geothermal_flux(%d, %d) = %f (%p)\n", m_grid->xs, m_grid->xs, rate.geothermal_flux(m_grid->xs, m_grid->xs), &rate.geothermal_flux(m_grid->xs, m_grid->xs));
rate.geothermal_flux.end_access();

cur.geothermal_flux.begin_access();
printf("GG cur.geothermal_flux(%d, %d) = %f (%p)\n", m_grid->xs, m_grid->xs, cur.geothermal_flux(m_grid->xs, m_grid->xs), &cur.geothermal_flux(m_grid->xs, m_grid->xs));
cur.geothermal_flux.end_access();

base.geothermal_flux.begin_access();
printf("GG base.geothermal_flux(%d, %d) = %f (%p)\n", m_grid->xs, m_grid->xs, base.geothermal_flux(m_grid->xs, m_grid->xs), &base.geothermal_flux(m_grid->xs, m_grid->xs));
base.geothermal_flux.end_access();
#endif
}

/** @param t0 Time of last time we coupled. */
void IBIceModel::prepare_outputs(double t0) {
  printf("BEGIN IBIceModel::prepare_outputs()\n");

  // ------ Difference between now and the last time we were called
  double t1 = enthalpy_t(); // Current time of the enthalpy portion of ice model.
  set_rate(t1 - t0);

  // ice_surface_enth & ice_surfac_enth_depth
  prepare_initial_outputs();

  printf("END IBIceModel::prepare_outputs()\n");
}

void IBIceModel::prepare_initial_outputs() {
  double ice_density = m_config->get_double("constants.ice.density", "kg m-3");

  const IceModelVec3 &ice_enthalpy = m_energy_model->enthalpy();

  AccessList access{ &ice_enthalpy, &M1, &M2, &H1, &H2, &V1, &V2, &m_geometry.ice_thickness };
  for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
    for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
      double const *Enth = ice_enthalpy.get_column(i, j);

      // Top Layer
      int const ks = m_grid->kBelowHeight(m_geometry.ice_thickness(i, j));
      V1(i, j) = m_geometry.ice_thickness(i, j) - m_grid->z(ks); // [m^3 m-2]
      M1(i, j) = V1(i, j) * ice_density;                // [kg m-2] = [m^3 m-2] [kg m-3]
      H1(i, j) = Enth[ks] * M1(i, j);                   // [J m-2] = [J kg-1] [kg m-2]

      // Second layer
      int const ks2 = ks - 1;
      if (ks2 >= 0) {
        V2(i, j) = m_grid->z(ks) - m_grid->z(ks2);
        M2(i, j) = V2(i, j) * ice_density;
        H2(i, j) = Enth[ks2] * M2(i, j);
      } else {
        // There is no second layer
        V2(i, j) = 0;
        M2(i, j) = 0;
        H2(i, j) = 0;
      }
    }
  }

  // ====================== Write to the post_energy.nc file (OPTIONAL)
}

void IBIceModel::time_setup() {
  // super::m_grid_setup() trashes m_time->start().  Now set it correctly.
  m_time->set_start(params.time_start_s);
  m_time->set(params.time_start_s);

  m_log->message(2, "* Run time: [%s, %s]  (%s years, using the '%s' calendar)\n", m_time->start_date().c_str(),
                 m_time->end_date().c_str(), m_time->run_length().c_str(), m_time->calendar().c_str());
}


void IBIceModel::misc_setup() {
  super::misc_setup();


  // ------ Initialize MassEnth structures: base, cur, rate
  for (auto &ii : cur.all_vecs) {
    ii.vec.set(0);
  }
  compute_enth2(cur.total.enth, cur.total.mass);
  cur.set_epsilon(m_grid);

  // base = cur
  auto base_ii(base.all_vecs.begin());
  auto cur_ii(cur.all_vecs.begin());
  for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii) {
    base_ii->vec.copy_from(cur_ii->vec);
  }
}

/** Sums over columns to compute enthalpy on 2D m_grid->

NOTE: Unfortunately so far PISM does not keep track of enthalpy in
"partially-filled" cells, so Enth2(i,j) is not valid at locations like
this one. We need to address this, but so far, it seems to me, the
best thing we can do is estimate Enth2(i,j) at partially-filled cells
by computing the average over icy neighbors. I think you can re-use
the idea from IceModel::get_threshold_thickness(...) (iMpartm_grid->cc).  */


void IBIceModel::compute_enth2(pism::IceModelVec2S &enth2, pism::IceModelVec2S &mass2) {
  //   getInternalColumn() is allocated already
  double ice_density = m_config->get_double("constants.ice.density", "kg m-3");

  const IceModelVec3 *ice_enthalpy = &m_energy_model->enthalpy();

  AccessList access{ &m_geometry.ice_thickness, ice_enthalpy, &enth2, &mass2 };
  for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
    for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
      enth2(i, j) = 0;
      mass2(i, j) = 0;

      // count all ice, including cells that have so little they
      // are considered "ice-free"
      if (m_geometry.ice_thickness(i, j) > 0) {
        const int ks       = m_grid->kBelowHeight(m_geometry.ice_thickness(i, j));
        double const *Enth = ice_enthalpy->get_column(i, j);
        for (int k = 0; k < ks; ++k) {
          double dz = (m_grid->z(k + 1) - m_grid->z(k));
          enth2(i, j) += Enth[k] * dz; // m J / kg
        }

        // Do the last layer a bit differently
        double dz = (m_geometry.ice_thickness(i, j) - m_grid->z(ks));
        enth2(i, j) += Enth[ks] * dz;
        enth2(i, j) *= ice_density;                        // --> J/m^2
        mass2(i, j) = m_geometry.ice_thickness(i, j) * ice_density; // --> kg/m^2
      }
    }
  }
}


/** Merges surface temperature derived from the energy balance model into any NaN values
in the vector provided.
@param deltah IN: Input from Icebin (change in enthalpy of each m_grid
    cell over the timestep) [W m-2].
@param default_val: The value that deltah(i,j) will have if no value
    is listed for that m_grid cell
@param timestep_s: Length of the current coupling timestep [s]
@param surface_temp OUT: Resulting surface temperature to use as the Dirichlet B.C.
*/
void IBIceModel::construct_surface_temp(
    pism::IceModelVec2S &deltah, // IN: Input from Icebin
    double default_val,
    double timestep_s,                 // Length of this coupling interval [s]
    pism::IceModelVec2S &surface_temp) // OUT: Temperature @ top of ice sheet (to use for Dirichlet B.C.)

{
  printf("BEGIN IBIceModel::merge_surface_temp default_val=%g\n", default_val);
  EnthalpyConverter::Ptr EC = ctx()->enthalpy_converter();

  double ice_density = m_config->get_double("constants.ice.density");

  const IceModelVec3 &ice_enthalpy = m_energy_model->enthalpy();

  {
    AccessList access{ &ice_enthalpy, &deltah, &m_geometry.ice_thickness, &surface_temp };

    // First time around, set effective_surface_temp to top temperature
    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        double &surface_temp_ij(surface_temp(i, j));
        double const &deltah_ij(deltah(i, j));

        double const *Enth = ice_enthalpy.get_column(i, j);

        // Enthalpy at top of ice sheet
        const int ks      = m_grid->kBelowHeight(m_geometry.ice_thickness(i, j));
        double spec_enth3 = Enth[ks]; // Specific enthalpy [J kg-1]

        if (deltah_ij != default_val) {
          // Adjust enthalpy @top by deltah
          double toplayer_dz = m_geometry.ice_thickness(i, j) - m_grid->z(ks); // [m]

          // [J kg-1] = [J kg-1]
          //     + [J m-2 s-1] * [m^2 m-3] * [m^3 kg-1] * [s]
          spec_enth3 = spec_enth3 + deltah_ij / (toplayer_dz * ice_density) * timestep_s;
        }


        // Convert specific enthalpy value to surface temperature
        const double p  = 0.0;                            // Pressure at top of ice sheet
        surface_temp_ij = EC->temperature(spec_enth3, p); // [K]
      }
    }
  }

  printf("END IBIceModel::merge_surface_temp\n");
}
}
} // namespace icebin::gpism
