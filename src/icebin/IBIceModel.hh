#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>

#include <pism/icemodel/IceModel.hh>
#include <pism/util/IceGrid.hh>

#include <pism/util/pism_options.hh>
#include <pism/coupler/atmosphere/Factory.hh>
#include <pism/coupler/ocean/Factory.hh>
#include <pism/coupler/surface/Factory.hh>
#include <pism/hydrology/NullTransport.hh>

#include <pism/util/Time.hh>
// --------------------------------
#include <pism/icebin/IBSurfaceModel.hh>
#include <pism/icebin/MassEnergyBudget.hh>

// Stuff defined in the icebin library
// (NOT a dependency of ours)
namespace icebin {
namespace gpism {
class IceModel_PISM;
}
}

namespace pism {
namespace icebin {


/** This is the IceBin customized version of PISM's pism::IceModel class.

See https://github.com/pism/pism/issues/219

Here's a short term solution, though: create a new class
IBIceModel derived from IceModel and re-implement
IceModel::allocate_couplers(). In it, set
IceModel::external_surface_model and IceModel::external_ocean_model as
you see fit (IceModel will not de-allocate a surface (ocean) model if
it is set to true) and allocate PSConstantICEBIN. You might also want
to add IBIceModel::get_surface_model() which returns
IceModel::surface to get access to PSConstantICEBIN from outside of
IBIceModel.
*/

class IBIceModel : public pism::IceModel {
  friend class ::icebin::gpism::IceModel_PISM;

public:
  typedef pism::IceModel super;
  struct Params {
    double time_start_s;
    std::string output_dir;
  };
  Params const params;

protected:
  MassEnergyBudget base; // Cumulative totals at start of this timestep
  MassEnergyBudget cur;  // Cumulative totals now
  MassEnergyBudget rate; // At end of coupling timestep, set to (cur - base) / dt

  // Output variables prepared for return to GCM
  // (relevant ice model state to be exported)
  pism::IceModelVec2S M1, M2;
  pism::IceModelVec2S H1, H2;
  pism::IceModelVec2S V1, V2;

protected:
  // see iceModel.cc
  virtual void allocate_storage();

public:
  virtual void massContExplicitStep(double dt,
                                    const IceModelVec2Stag &diffusive_flux,
                                    const IceModelVec2V &advective_velocity);
  virtual void accumulateFluxes_massContExplicitStep(int i, int j,
                                                     double surface_mass_balance, // [m s-1] ice equivalent (from PISM)
                                                     double meltrate_grounded,    // [m s-1] ice equivalent
                                                     double meltrate_floating,    // [m s-1] ice equivalent
                                                     double divQ_SIA,             // [m s-1] ice equivalent
                                                     double divQ_SSA,             // [m s-1] ice equivalent
                                                     double Href_to_H_flux,       // [m s-1] ice equivalent
                                                     double nonneg_rule_flux);    // [m s-1] ice equivalent
private:
  // Temporary variables inside massContExplicitStep()
  double _ice_density;              // From config
  double _meter_per_s_to_kg_per_m2; // Conversion factor computed from _ice_density


private:
  // Utility function
  void prepare_nc(std::string const &fname, std::unique_ptr<pism::PIO> &nc);

public:
  /** @param t0 Time of last time we coupled. */
  void set_rate(double dt);

  void reset_rate();


  std::unique_ptr<pism::PIO> pre_mass_nc; //!< Write variables every time massContPostHook() is called.
  std::unique_ptr<pism::PIO> post_mass_nc;
  std::unique_ptr<pism::PIO> pre_energy_nc;
  std::unique_ptr<pism::PIO> post_energy_nc;

  // see iceModel.cc for implementation of constructor and destructor:
  /** @param gcm_params Pointer to IceModel::gcm_params.  Lives at least as long as this object. */
  IBIceModel(IceGrid::Ptr g, Context::Ptr context, IBIceModel::Params const &_params);
  virtual ~IBIceModel(); // must be virtual merely because some members are virtual

  virtual void allocate_subglacial_hydrology();
  virtual void allocate_couplers();
  virtual void time_setup();
  virtual void misc_setup();

  void compute_enth2(pism::IceModelVec2S &enth2, pism::IceModelVec2S &mass2);

  /** @return Our instance of IBSurfaceModel */
  pism::icebin::IBSurfaceModel *ib_surface_model() {
    return dynamic_cast<IBSurfaceModel *>(m_surface.get());
  }
  pism::hydrology::NullTransport* null_hydrology() {
    return dynamic_cast<hydrology::NullTransport *>(pism::IceModel::m_subglacial_hydrology.get());
  }


  /** @return Current time for mass timestepping */
  double mass_t() const {
    return m_time->current();
  }
  /** @return Current time for enthalpy timestepping */
  double enthalpy_t() const {
    return t_TempAge;
  }

  void energy_step();

  void prepare_outputs(double time_s);

  /** Read things out of the ice model that will be sent back BEFORE
    the first coupling timestep (eg, ice surface enthalpy) */
  void prepare_initial_outputs();

  /** Merges surface temperature derived from the energy balance model into any NaN values
    in the vector provided.
    @param deltah IN: Input from Icebin (change in enthalpy of each grid
        cell over the timestep) [W m-2].
    @param default_val: The value that deltah(i,j) will have if no value
        is listed for that grid cell
    @param timestep_s: Length of the current coupling timestep [s]
    @param surface_temp OUT: Resulting surface temperature to use as the Dirichlet B.C.
    */
  void construct_surface_temp(pism::IceModelVec2S &deltah, // IN: Input from Icebin
                              double default_val,
                              double timestep_s, // Length of this coupling interval [s]
                              pism::IceModelVec2S &surface_temp);
};
}
}
