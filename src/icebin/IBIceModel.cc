#include <iostream>

#include <base/energy/bedrockThermalUnit.hh>
#include <base/enthalpyConverter.hh>
#include <base/util/io/PIO.hh>
#include <base/util/io/io_helpers.hh>

#include <icebin/IBIceModel.hh>
#include <icebin/IBSurfaceModel.hh>


namespace pism {
namespace icebin {

// ================================

IBIceModel::IBIceModel(IceGrid::Ptr g, Context::Ptr context, IBIceModel::Params const &_params) :
    pism::IceModel(g, context),
    params(_params) {
  // empty
}

IBIceModel::~IBIceModel() {
  // empty
}


void IBIceModel::allocate_subglacial_hydrology() {
    printf("BEGIN IBIceModel::allocate_subglacial_hydrology()\n");
    if (pism::IceModel::subglacial_hydrology) return; // indicates it has already been allocated
    subglacial_hydrology = new pism::icebin::NullTransportHydrology(m_grid);
    printf("END IBIceModel::allocate_subglacial_hydrology()\n");
}


void IBIceModel::allocate_couplers() {
  // Initialize boundary models:
  atmosphere::Factory pa(m_grid);
  surface::Factory ps(m_grid);
  ocean::Factory po(m_grid);
  atmosphere::AtmosphereModel *atmosphere;

  if (m_surface == NULL) {

    m_log->message(2,
             "# Allocating a surface process model or coupler...\n");

    m_surface = new IBSurfaceModel(m_grid);
    m_external_surface_model = false;

    atmosphere = pa.create();
    m_surface->attach_atmosphere_model(atmosphere);
  }

  if (m_ocean == NULL) {
    m_log->message(2,
             "# Allocating an ocean model or coupler...\n");

    m_ocean = po.create();
    m_external_ocean_model = false;
  }
}


void IBIceModel::createVecs() {
    super::createVecs();

    printf("BEGIN IBIceModel::createVecs()\n");
    base.create(m_grid, "", WITHOUT_GHOSTS);
    cur.create(m_grid, "", WITHOUT_GHOSTS);
    rate.create(m_grid, "", WITHOUT_GHOSTS);
    printf("END IBIceModel::createVecs()\n");

    M1.create(m_grid, "M1", pism::WITHOUT_GHOSTS);
    M2.create(m_grid, "M2", pism::WITHOUT_GHOSTS);
    M1.create(m_grid, "H1", pism::WITHOUT_GHOSTS);
    M2.create(m_grid, "H2", pism::WITHOUT_GHOSTS);
    M1.create(m_grid, "V1", pism::WITHOUT_GHOSTS);
    M2.create(m_grid, "V2", pism::WITHOUT_GHOSTS);

    std::cout << "IBIceModel Conservation Formulas:" << std::endl;
    cur.print_formulas(std::cout);
}

void IBIceModel::massContPreHook() {
#if 0   // Avoid warnings until we put something in this method

    // Enthalpy and mass continuity are stepped with different timesteps.
    // Fish out the timestep relevant to US.
    const double my_t0 = m_grid.time->current();
    const double my_dt = this->dt;
#endif
}


void IBIceModel::massContPostHook() {
#if 0   // Avoid warnings until we put something in this method

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
        const IceModelVec2S &upward_geothermal_flux(btu->upward_geothermal_flux());
        cur.upward_geothermal_flux.add(my_dt, upward_geothermal_flux);
    }

    // ----------- Geothermal Flux
    cur.geothermal_flux.add(my_dt, m_geothermal_flux);

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

    _ice_density = m_config->get_double("ice_density");
    _meter_per_s_to_kg_per_m2 = m_dt * _ice_density;


    // =========== The Mass Continuity Step Itself
    // This will call through to accumulateFluxes_massContExplicitStep()
    // in the inner loop
    {
      AccessList access{&cur.pism_smb, &cur.melt_grounded,
          &cur.melt_floating, &cur.internal_advection, &cur.href_to_h,
          &cur.nonneg_rule};

      super::massContExplicitStep();
    }

    // =========== AFTER the Mass Continuity Step

    // ----------- SMB: Pass inputs through to outputs.
    // They are needed to participate in mass/energy budget
    IBSurfaceModel *ib_surface = ib_surface_model();


    {
      AccessList access{&ib_surface->icebin_massxfer,
          &ib_surface->icebin_enthxfer,
          &ib_surface->icebin_deltah,
          &cur.icebin_xfer,
          &cur.icebin_deltah};

      for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
        for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
          cur.icebin_xfer.mass(i,j) += m_dt * ib_surface->icebin_massxfer(i,j);
          cur.icebin_xfer.enth(i,j) += m_dt * ib_surface->icebin_enthxfer(i,j);
          cur.icebin_deltah(i,j) += m_dt * ib_surface->icebin_deltah(i,j);
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
    double surface_mass_balance,        // [m s-1] ice equivalent
    double meltrate_grounded,           // [m s-1] ice equivalent
    double meltrate_floating,           // [m s-1] ice equivalent
    double divQ_SIA,                    // [m s-1] ice equivalent
    double divQ_SSA,                    // [m s-1] ice equivalent
    double Href_to_H_flux,              // [m] ice equivalent
    double nonneg_rule_flux)            // [m] ice equivalent
{
    EnthalpyConverter::Ptr EC = ctx()->enthalpy_converter();

    // -------------- Melting
    double p_basal = EC->pressure(m_ice_thickness(i,j));
    double T = EC->melting_temperature(p_basal);
    double specific_enth_basal = EC->enthalpy_permissive(T, 1.0, p_basal);
    double mass;

    // ------- Melting at base of ice sheet
    mass = -meltrate_grounded * _meter_per_s_to_kg_per_m2;
    cur.melt_grounded.mass(i,j) += mass;
    cur.melt_grounded.enth(i,j) += mass * specific_enth_basal;

    // ------- Melting under ice shelf
    mass = -meltrate_floating * _meter_per_s_to_kg_per_m2;
    cur.melt_floating.mass(i,j) += mass;
    cur.melt_floating.enth(i,j) += mass * specific_enth_basal;


    // -------------- internal_advection
    const int ks = m_grid->kBelowHeight(m_ice_thickness(i,j));
    double *Enth = m_ice_enthalpy.get_column(i,j);
    double specific_enth_top = Enth[ks];        // Approximate, we will use the enthalpy of the top layer...

    mass = -(divQ_SIA + divQ_SSA) * _meter_per_s_to_kg_per_m2;

    cur.internal_advection.mass(i,j) += mass;
    cur.internal_advection.enth(i,j) += mass * specific_enth_top;


    // -------------- Get the easy veriables out of the way...
    mass = surface_mass_balance * _meter_per_s_to_kg_per_m2;
    cur.pism_smb.mass(i,j) += mass;
    cur.pism_smb.enth(i,j) += mass * specific_enth_top;
    cur.nonneg_rule(i,j) -= nonneg_rule_flux * _ice_density;
    cur.href_to_h(i,j) += Href_to_H_flux * _ice_density;


//  printf("END IBIceModel::accumulateFluxes_MassContExplicitStep()\n");
}



void IBIceModel::prepare_nc(std::string const &fname, std::unique_ptr<PIO> &nc) {

//    nc.reset(new PIO(m_grid->com, m_grid->ctx()->config()->get_string("output_format")));

    nc.reset(new PIO(m_grid->com, m_config->get_string("output_format")));

    nc->open(fname, PISM_READWRITE_MOVE);
    io::define_time(*nc, m_grid->ctx()->config()->get_string("time_dimension_name"),
        m_grid->ctx()->time()->calendar(),
        m_grid->ctx()->time()->CF_units_string(),
        m_grid->ctx()->unit_system());

// These are in iMtimseries, but not listed as required in iceModelVec.hh
//    nc->put_att_text(m_config.get_string("time_dimension_name"),
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
        AccessList access{&vbase, &vcur, &vrate};
        for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
          for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
            // rate = cur - base: Just for DELTA and EPISLON flagged vectors
            if (base_ii->flags & (MassEnergyBudget::DELTA | MassEnergyBudget::EPSILON)) {
              vrate(i,j) = (vcur(i,j) - vbase(i,j)) * by_dt;
            } else {
              // Or else just copy the to rate
              vrate(i,j) = vcur(i,j);
            }
          }
        }
      }
    }

    printf("END IBIceModel::set_rate()\n");
}

void IBIceModel::reset_rate()
{
    // Compute differences, and set base = cur
    auto base_ii(base.all_vecs.begin());
    auto cur_ii(cur.all_vecs.begin());
    for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii) {
        IceModelVec2S &vbase(base_ii->vec);
        IceModelVec2S &vcur(cur_ii->vec);

        // This cannot go in the loop above with PETSc because
        // vbase is needed on the RHS of the equations above.
        AccessList access{&vbase, &vcur};
        for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
        for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
            // base = cur: For ALL vectors
            vbase(i,j) = vcur(i,j);
        }}
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
void IBIceModel::prepare_outputs(double t0)
{
    printf("BEGIN IBIceModel::prepare_outputs()\n");

    // ------ Difference between now and the last time we were called
    double t1 = enthalpy_t();   // Current time of the enthalpy portion of ice model.
    set_rate(t1 - t0);

    // ice_surface_enth & ice_surfac_enth_depth
    prepare_initial_outputs();

    // ------ Write it out
#if 0
    // This is not really needed, since Icebin also writes out
    // the same fields.
    PIO nc(m_grid, m_grid->m_config.get_string("output_format"));
    nc.open((params.output_dir / "post_energy.nc").c_str(), PISM_READWRITE);    // append to file
    nc.append_time(m_config.get_string("time_dimension_name"), t1);
    m_ice_enthalpy.write(nc, PISM_DOUBLE);
    ice_thickness.write(nc, PISM_DOUBLE);
    ice_surface_temp.write(nc, PISM_DOUBLE);
    PSConstantICEBIN *surface = ps_constant_icebin();
    surface->effective_surface_temp.write(nc, PISM_DOUBLE);
    for (auto ii = rate.all_vecs.begin(); ii != rate.all_vecs.end(); ++ii) {
        ii->vec.write(nc, PISM_DOUBLE);
    }
    nc.close();
#endif

    printf("END IBIceModel::prepare_outputs()\n");
}

void IBIceModel::prepare_initial_outputs()
{
    double ice_density = m_config->get_double("ice_density"); // [kg m-3]

    // --------- ice_surface_enth from m_ice_enthalpy
    AccessList access{&m_ice_enthalpy, &M1, &M2, &H1, &H2, &V1, &V2, &m_ice_thickness};
    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        double const *Enth = m_ice_enthalpy.get_column(i,j);

        // Top Layer
        int const ks = m_grid->kBelowHeight(m_ice_thickness(i,j));
        V1(i,j) = m_ice_thickness(i,j) - m_grid->z(ks);    // [m^3 m-2]
        M1(i,j) = V1(i,j) * ice_density;    // [kg m-2] = [m^3 m-2] [kg m-3]
        H1(i,j) = Enth[ks] * M1(i,j);       // [J m-2] = [J kg-1] [kg m-2]

        // Second layer
        int const ks2 = ks - 1;
        if (ks2 >= 0) {
          V2(i,j) = m_grid->z(ks) - m_grid->z(ks2);
          M2(i,j) = V2(i,j) * ice_density;
          H2(i,j) = Enth[ks2] * M2(i,j);
        } else {
          // There is no second layer
          V2(i,j) = 0;
          M2(i,j) = 0;
          H2(i,j) = 0;
        }
      }
    }

    // ====================== Write to the post_energy.nc file (OPTIONAL)
}

void IBIceModel::time_setup()
{
  // super::m_grid_setup() trashes m_time->start().  Now set it correctly.
  m_time->set_start(params.time_start_s);
  m_time->set(params.time_start_s);

  m_log->message(2,
             "* Run time: [%s, %s]  (%s years, using the '%s' calendar)\n",
             m_time->start_date().c_str(),
             m_time->end_date().c_str(),
             m_time->run_length().c_str(),
             m_time->calendar().c_str());
}


void IBIceModel::misc_setup()
{
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

#if 0
    // ---------- Create the netCDF output file
    std::unique_ptr<PIO> nc;
    std::string ofname = (params.output_dir / "post_energy.nc").string();
    prepare_nc(ofname, nc);

    // -------- Define MethEnth structres in netCDF file
    m_ice_enthalpy.define(*nc, PISM_DOUBLE);
    ice_thickness.define(*nc, PISM_DOUBLE);
    ice_surface_temp.define(*nc, PISM_DOUBLE);
    PSConstantICEBIN *surface = ps_constant_icebin();
    surface->effective_surface_temp.define(*nc, PISM_DOUBLE);
    for (auto ii = rate.all_vecs.begin(); ii != rate.all_vecs.end(); ++ii) {
        ii->vec.define(*nc, PISM_DOUBLE);
    }

    // --------- Close and return
    nc->close();
#endif
}

/** Sums over columns to compute enthalpy on 2D m_grid->

NOTE: Unfortunately so far PISM does not keep track of enthalpy in
"partially-filled" cells, so Enth2(i,j) is not valid at locations like
this one. We need to address this, but so far, it seems to me, the
best thing we can do is estimate Enth2(i,j) at partially-filled cells
by computing the average over icy neighbors. I think you can re-use
the idea from IceModel::get_threshold_thickness(...) (iMpartm_grid->cc).  */


void IBIceModel::compute_enth2(pism::IceModelVec2S &enth2, pism::IceModelVec2S &mass2)
{
    //   getInternalColumn() is allocated already
    double ice_density = m_config->get_double("ice_density");
    AccessList access{&m_ice_thickness, &m_ice_enthalpy, &enth2, &mass2};
    for (int i=m_grid->xs(); i<m_grid->xs() +m_grid->xm(); ++i) {
        for (int j=m_grid->ys(); j<m_grid->ys() + m_grid->ym(); ++j) {
            enth2(i,j) = 0;
            mass2(i,j) = 0;

            // count all ice, including cells that have so little they
            // are considered "ice-free"
            if (m_ice_thickness(i,j) > 0) {
                const int ks = m_grid->kBelowHeight(m_ice_thickness(i,j));
                double const *Enth = m_ice_enthalpy.get_column(i,j); // do NOT delete this pointer: space returned by
                for (int k=0; k<ks; ++k) {
                    double dz = (m_grid->z(k+1) - m_grid->z(k));
                    enth2(i,j) += Enth[k] * dz;     // m J / kg
                }

                // Do the last layer a bit differently
                double dz = (m_ice_thickness(i,j) - m_grid->z(ks));
                enth2(i,j) += Enth[ks] * dz;
                enth2(i,j) *= ice_density;      // --> J/m^2
                mass2(i,j) = m_ice_thickness(i,j) * ice_density;      // --> kg/m^2
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
    pism::IceModelVec2S &deltah,            // IN: Input from Icebin
    double default_val,
    double timestep_s,      // Length of this coupling interval [s]
    pism::IceModelVec2S &surface_temp)  // OUT: Temperature @ top of ice sheet (to use for Dirichlet B.C.)

{
  printf("BEGIN IBIceModel::merge_surface_temp default_val=%g\n", default_val);
  EnthalpyConverter::Ptr EC = ctx()->enthalpy_converter();

  double ice_density = m_config->get_double("ice_density");

  {
    AccessList access{&m_ice_enthalpy, &deltah, &m_ice_thickness, &surface_temp};

    // First time around, set effective_surface_temp to top temperature
    for (int i = m_grid->xs(); i < m_grid->xs() + m_grid->xm(); ++i) {
      for (int j = m_grid->ys(); j < m_grid->ys() + m_grid->ym(); ++j) {
        double &surface_temp_ij(surface_temp(i,j));
        double const &deltah_ij(deltah(i,j));

        double const *Enth = m_ice_enthalpy.get_column(i,j);

        // Enthalpy at top of ice sheet
        const int ks = m_grid->kBelowHeight(m_ice_thickness(i,j));
        double spec_enth3 = Enth[ks];       // Specific enthalpy [J kg-1]

        if (deltah_ij != default_val) {
          // Adjust enthalpy @top by deltah
          double toplayer_dz = m_ice_thickness(i,j) - m_grid->z(ks);     // [m]

          // [J kg-1] = [J kg-1]
          //     + [J m-2 s-1] * [m^2 m-3] * [m^3 kg-1] * [s]
          spec_enth3 = spec_enth3
            + deltah_ij / (toplayer_dz * ice_density) * timestep_s;
        }


        // Convert specific enthalpy value to surface temperature
        const double p = 0.0;       // Pressure at top of ice sheet
        surface_temp_ij = EC->temperature(spec_enth3, p);    // [K]
      }
    }
  }

  printf("END IBIceModel::merge_surface_temp\n");
}


}}  // namespace icebin::gpism
