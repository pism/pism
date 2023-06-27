#include <iostream>
#include <limits>

#include <base/energy/BedThermalUnit.hh>
#include <base/enthalpyConverter.hh>
#include <base/util/io/PIO.hh>
#include <base/util/io/io_helpers.hh>
#include <base/util/Mask.hh>
#include <earth/PISMBedDef.hh>

#include <icebin/IBIceModel.hh>
#include <icebin/IBSurfaceModel.hh>


namespace pism {
namespace icebin {

static double const NaN = std::numeric_limits<double>::quiet_NaN();

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
  if (pism::IceModel::m_subglacial_hydrology)
    return; // indicates it has already been allocated
  m_subglacial_hydrology = new pism::icebin::NullTransportHydrology(m_grid);
  printf("END IBIceModel::allocate_subglacial_hydrology()\n");
}


void IBIceModel::allocate_couplers() {
  // Initialize boundary models:
  atmosphere::Factory pa(m_grid);
  surface::Factory ps(m_grid);
  ocean::Factory po(m_grid);
  atmosphere::AtmosphereModel *atmosphere;

  if (m_surface == NULL) {

    m_log->message(2, "# Allocating a surface process model or coupler...\n");

    m_surface                = new IBSurfaceModel(m_grid);

    atmosphere = pa.create();
    m_surface->attach_atmosphere_model(atmosphere);
  }

  if (m_ocean == NULL) {
    m_log->message(2, "# Allocating an ocean model or coupler...\n");

    m_ocean                = po.create();
  }
}


void IBIceModel::createVecs() {
  super::createVecs();

  printf("BEGIN IBIceModel::createVecs()\n");
  base.create(m_grid, "", WITHOUT_GHOSTS);
  cur.create(m_grid, "", WITHOUT_GHOSTS);
  rate.create(m_grid, "", WITHOUT_GHOSTS);
  printf("END IBIceModel::createVecs()\n");

  // Sent back to IceBin
  ice_top_senth.create(m_grid, "ice_top_senth", pism::WITHOUT_GHOSTS);
  elevmask_ice.create(m_grid, "elevmask_ice", pism::WITHOUT_GHOSTS);
  elevmask_land.create(m_grid, "elevmask_land", pism::WITHOUT_GHOSTS);

  std::cout << "IBIceModel Conservation Formulas:" << std::endl;
  cur.print_formulas(std::cout);
}

void IBIceModel::massContPreHook() {
#if 0 // Avoid warnings until we put something in this method

    // Enthalpy and mass continuity are stepped with different timesteps.
    // Fish out the timestep relevant to US.
    const double my_t0 = m_grid.time->current();
    const double my_dt = this->dt;
#endif
}


void IBIceModel::massContPostHook() {
#if 0 // Avoid warnings until we put something in this method

    // Enthalpy and mass continuity are stepped with different timesteps.
    // Fish out the timestep relevant to US.
    const double my_t0 = m_grid.time->current();
    const double my_dt = this->dt;
#endif
}


void IBIceModel::energyStep() {

  printf("BEGIN IBIceModel::energyStep(t=%f, dt=%f)\n", t_TempAge, dt_TempAge);

  // Enthalpy and mass continuity are stepped with different timesteps.
  // Fish out the timestep relevant to US.
  // const double my_t0 = t_TempAge;          // Time at beginning of timestep
  const double my_dt = dt_TempAge;

  // =========== BEFORE Energy Step

  // =========== The Energy Step Itself
  super::energyStep();

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

  printf("END IBIceModel::energyStep(time=%f)\n", t_TempAge);
}

void IBIceModel::massContExplicitStep() {

  printf("BEGIN IBIceModel::MassContExplicitStep()\n");

  _ice_density              = m_config->get_double("constants.ice.density");
  _meter_per_s_to_kg_per_m2 = m_dt * _ice_density;


  // =========== The Mass Continuity Step Itself
  // This will call through to accumulateFluxes_massContExplicitStep()
  // in the inner loop.  We must open access to variables we will use
  // in that subroutine.
  {
    AccessList access{ &cur.pism_smb,           &cur.melt_grounded, &cur.melt_floating,
                       &cur.internal_advection, &cur.href_to_h,     &cur.nonneg_rule };

    super::massContExplicitStep();
  }

  // =========== AFTER the Mass Continuity Step

  // ----------- SMB: Pass inputs through to outputs.
  // They are needed to participate in mass/energy budget
  // This allows them to be used in the MassEnergyBudget computation.
  IBSurfaceModel *ib_surface = ib_surface_model();

  {
    AccessList access{ &ib_surface->massxfer, &ib_surface->enthxfer, &ib_surface->deltah,
                       &cur.smb, &cur.deltah };

    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        cur.smb.mass(i, j) += m_dt * ib_surface->massxfer(i, j);
        cur.smb.enth(i, j) += m_dt * ib_surface->enthxfer(i, j);
        cur.deltah(i, j) += m_dt * ib_surface->deltah(i, j);
      }
    }
  }


  printf("END IBIceModel::MassContExplicitStep()\n");
}


/** This is called IMMEDIATELY after ice is gained/lost in
iMgeometry.cc (massContExplicitStep()).  Here we can record the same
values that PISM saw when moving ice around. */
void IBIceModel::accumulateFluxes_massContExplicitStep(
  int i, int j,
  double surface_mass_balance, // [m s-1] ice equivalent
  double basal_melt_rate,      // [m s-1] ice equivalent
  double divQ_SIA,             // [m s-1] ice equivalent
  double divQ_SSA,             // [m s-1] ice equivalent
  double Href_to_H_flux,       // [m] ice equivalent
  double nonneg_rule_flux)     // [m] ice equivalent
{
// printf("BEGIN IBIceModel::accumulateFluxes_MassContExplicitStep()\n");

  EnthalpyConverter::Ptr EC = ctx()->enthalpy_converter();

  // -------------- Melting
  double p_basal             = EC->pressure(m_ice_thickness(i, j));
  double T                   = EC->melting_temperature(p_basal);
  double specific_enth_basal = EC->enthalpy_permissive(T, 1.0, p_basal);
  double mass;

  // ------- Melting at base of ice sheet
  mass = -basal_melt_rate * _meter_per_s_to_kg_per_m2;
  // TODO: Change this to just cur.melt_rate
  cur.melt_grounded.mass(i, j) += mass;
  cur.melt_grounded.enth(i, j) += mass * specific_enth_basal;

  // -------------- internal_advection
  const int ks             = m_grid->kBelowHeight(m_ice_thickness(i, j));
  double *Enth             = m_ice_enthalpy.get_column(i, j);
  double specific_enth_top = Enth[ks]; // Approximate, we will use the enthalpy of the top layer...

  mass = -(divQ_SIA + divQ_SSA) * _meter_per_s_to_kg_per_m2;

  cur.internal_advection.mass(i, j) += mass;
  cur.internal_advection.enth(i, j) += mass * specific_enth_top;


  // -------------- Get the easy variables out of the way...
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

/** Differences and divides by accumulated time over coupling timestep.
Eg, to convert [kg m-2] --> [kg m-2 s-1]
@param t0 Time of last time we coupled. */
void IBIceModel::set_rate(double dt) {

  printf("BEGIN IBIceModel::set_rate(dt=%f)\n", dt);

  if (dt == 0) throw RuntimeError(PISM_ERROR_LOCATION,
    "Coupling timestep has size dt=0");

  double by_dt = 1.0 / dt;
printf("IBIceModel::set_rate(dt=%f)\n", dt);

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
            // Or else just copy the to "rate"
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
}

void IBIceModel::prepare_outputs(double time_s) {
  EnthalpyConverter::Ptr EC = ctx()->enthalpy_converter();

  // --------- ice_surface_enth from m_ice_enthalpy
  auto &ice_surface_elevation(this->ice_surface_elevation());
  auto &bed_topography(this->bed_model()->bed_elevation());

  auto &cell_type(this->cell_type());
  AccessList access{
    &m_ice_enthalpy, &m_ice_thickness,        // INPUTS
    &ice_surface_elevation, &cell_type,
    &ice_top_senth, &elevmask_ice, &elevmask_land };                 // OUTPUT
  for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
    for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
      double const *Enth = m_ice_enthalpy.get_column(i, j);

      // Top Layer
      int const ks = m_grid->kBelowHeight(m_ice_thickness(i, j));
      double senth = Enth[ks];   // [J kg-1]
      ice_top_senth(i,j) = senth;

      // elevmask_ice and elevmask_land: Used by IceBin for elevation and masking
      switch((int)cell_type(i,j)) {
        case MASK_GROUNDED:
        case MASK_FLOATING:
          elevmask_ice(i,j) = ice_surface_elevation(i,j);
          elevmask_land(i,j) = ice_surface_elevation(i,j);
        break;
        case MASK_ICE_FREE_BEDROCK:
          elevmask_ice(i,j) = NaN;
          elevmask_land(i,j) = ice_surface_elevation(i,j);
        break;
        case MASK_ICE_FREE_OCEAN:
        case MASK_UNKNOWN :
          elevmask_ice(i,j) = NaN;
          elevmask_land(i,j) = NaN;
        break;
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

/** Sums over columns to compute enthalpy on 2D m_grid.
This is used to track conservation of energy.

NOTE: Unfortunately so far PISM does not keep track of enthalpy in
"partially-filled" cells, so Enth2(i,j) is not valid at locations like
this one. We need to address this, but so far, it seems to me, the
best thing we can do is estimate Enth2(i,j) at partially-filled cells
by computing the average over icy neighbors. I think you can re-use
the idea from IceModel::get_threshold_thickness(...) (iMpartm_grid->cc).  */

void IBIceModel::compute_enth2(pism::IceModelVec2S &enth2, pism::IceModelVec2S &mass2) {
  //   getInternalColumn() is allocated already
  double ice_density = m_config->get_double("constants.ice.density", "kg m-3");
  AccessList access{
    &m_ice_thickness, &m_ice_enthalpy,    // Inputs
    &enth2, &mass2 };                     // Outputs
  for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
    for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
      enth2(i, j) = 0;
      mass2(i, j) = 0;

      // count all ice, including cells that have so little they
      // are considered "ice-free"
      if (m_ice_thickness(i, j) > 0) {
        const int ks       = m_grid->kBelowHeight(m_ice_thickness(i, j));
        double const *Enth = m_ice_enthalpy.get_column(i, j); // do NOT delete this pointer: space returned by
        for (int k = 0; k < ks; ++k) {
          double dz = (m_grid->z(k + 1) - m_grid->z(k));
          enth2(i, j) += Enth[k] * dz; // [m J kg-1]
        }

        // Last layer is (potentially) partial
        double dz = (m_ice_thickness(i, j) - m_grid->z(ks));
        enth2(i, j) += Enth[ks] * dz;    // [m J kg-1]

        // Finish off after all layers processed
        enth2(i, j) *= ice_density;                        // --> [J m-2]
        mass2(i, j) = m_ice_thickness(i, j) * ice_density; // --> [kg m-2]
      }
    }
  }
}


/** Merges surface temperature derived from m_ice_enthalpy into any NaN values
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

  {
    AccessList access{ &m_ice_enthalpy, &deltah, &m_ice_thickness, &surface_temp };

    // First time around, set effective_surface_temp to top temperature
    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        double &surface_temp_ij(surface_temp(i, j));
        double const &deltah_ij(deltah(i, j));

        double const *Enth = m_ice_enthalpy.get_column(i, j);

        // Enthalpy at top of ice sheet
        const int ks      = m_grid->kBelowHeight(m_ice_thickness(i, j));
        double spec_enth3 = Enth[ks]; // Specific enthalpy [J kg-1]

        if (deltah_ij != default_val) {
          // Adjust enthalpy @top by deltah
          double toplayer_dz = m_ice_thickness(i, j) - m_grid->z(ks); // [m]

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
