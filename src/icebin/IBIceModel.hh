#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>

#include <base/iceModel.hh>
#include <base/util/IceGrid.hh>

#include <base/util/pism_options.hh>
#include <coupler/atmosphere/PAFactory.hh>
#include <coupler/ocean/POFactory.hh>
#include <coupler/surface/PSFactory.hh>

#include <base/util/PISMTime.hh>
// --------------------------------
#include <icebin/IBSurfaceModel.hh>
#include <icebin/MassEnergyBudget.hh>
#include <icebin/NullTransportHydrology.hh>

// Stuff defined in the icebin library
// (NOT a dependency of ours)
namespace icebin {
namespace gpism {
class IceCoupler_PISM;
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
  friend class ::icebin::gpism::IceCoupler_PISM;

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

  // Specific enthalpy at surface of the ice sheet [J kg-1]
  pism::IceModelVec2S ice_top_senth;

public:
  // Elevation of ice grid cells, with NaN off the ice sheet [m]
  pism::IceModelVec2S elevmask_ice;
  // Elevation of ice+bare land grid cells, with NaN in the ocean [m]
  pism::IceModelVec2S elevmask_land;

protected:
  // see iceModel.cc
  virtual void createVecs();

public:
  virtual void accumulateFluxes_massContExplicitStep(
    int i, int j,
    double surface_mass_balance, // [m s-1] ice equivalent (from PISM)
    double basal_melt_rate,      // [m s-1] ice equivalent
    double divQ_SIA,             // [m s-1] ice equivalent
    double divQ_SSA,             // [m s-1] ice equivalent
    double Href_to_H_flux,       // [m s-1] ice equivalent
    double nonneg_rule_flux);    // [m s-1] ice equivalent
  virtual void massContExplicitStep();

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
    return dynamic_cast<IBSurfaceModel *>(m_surface);
  }
  pism::icebin::NullTransportHydrology *null_hydrology() {
    return dynamic_cast<NullTransportHydrology *>(pism::IceModel::m_subglacial_hydrology);
  }


  /** @return Current time for mass timestepping */
  double mass_t() const {
    return m_time->current();
  }
  /** @return Current time for enthalpy timestepping */
  double enthalpy_t() const {
    return t_TempAge;
  }

  // I added these...
  void massContPreHook();
  void massContPostHook();
  // Pre and post for energy
  void energyStep();

  void prepare_outputs(double time_s);

  /** Read things out of the ice model that will be sent back BEFORE
    the first coupling timestep (eg, ice surface enthalpy) */
  void prepare_initial_outputs();

  /** Merges surface temperature derived from m_ice_enthalpy into any NaN values
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
