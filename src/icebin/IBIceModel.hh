#pragma once

#include "pism/icebin/IBSurfaceModel.hh"
#include "pism/icebin/MassEnergyBudget.hh"
#include "pism/icebin/NullTransportHydrology.hh"
#include "pism/icemodel/IceModel.hh"
#include <memory>

// Stuff defined in the icebin library
// (NOT a dependency of ours)
namespace icebin {
namespace gpism {
class IceCoupler_PISM;
}
} // namespace icebin

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
  MassEnergyBudget base; // Cumulative totals at start of this time step
  MassEnergyBudget cur;  // Cumulative totals now
  MassEnergyBudget rate; // At end of coupling timestep, set to (cur - base) / dt

  // Output variables prepared for return to GCM
  // (relevant ice model state to be exported)

  // Specific enthalpy at surface of the ice sheet [J kg-1]
  pism::array::Scalar ice_top_senth;

public:
  // Elevation of ice grid cells, with NaN off the ice sheet [m]
  pism::array::Scalar elevmask_ice;
  // Elevation of ice+bare land grid cells, with NaN in the ocean [m]
  pism::array::Scalar elevmask_land;

protected:

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
  double m_ice_density;              // From config
  double m_meter_per_s_to_kg_per_m2; // Conversion factor computed from m_ice_density

public:
  /** @param t0 Time of last time we coupled. */
  void set_rate(double dt);

  void reset_rate();

  std::unique_ptr<pism::File> pre_mass_nc; //!< Write variables every time massContPostHook() is called.
  std::unique_ptr<pism::File> post_mass_nc;
  std::unique_ptr<pism::File> pre_energy_nc;
  std::unique_ptr<pism::File> post_energy_nc;

  // see iceModel.cc for implementation of constructor and destructor:
  /** @param gcm_params Pointer to IceModel::gcm_params.  Lives at least as long as this object. */
  IBIceModel(std::shared_ptr<pism::Grid> grid, const std::shared_ptr<Context> &context,
             IBIceModel::Params const &_params);
  virtual ~IBIceModel(); // must be virtual merely because some members are virtual

  virtual void allocate_subglacial_hydrology();
  virtual void allocate_couplers();
  virtual void time_setup();
  virtual void misc_setup();

  void compute_enth2(pism::array::Scalar &enth2, pism::array::Scalar &mass2);

  /** @return Our instance of IBSurfaceModel */
  pism::icebin::IBSurfaceModel *ib_surface_model() {
    return std::dynamic_pointer_cast<IBSurfaceModel>(m_surface).get();
  }

  pism::icebin::NullTransportHydrology *null_hydrology() {
    return dynamic_cast<pism::icebin::NullTransportHydrology *>(m_subglacial_hydrology.get());
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
  void energy_step();

  void prepare_outputs(double time_s);

  void dumpToFile(const std::string &filename) const ;

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
  void construct_surface_temp(pism::array::Scalar &deltah, // IN: Input from Icebin
                              double default_val,
                              double timestep_s, // Length of this coupling interval [s]
                              pism::array::Scalar &surface_temp);
};
} // namespace icebin
} // namespace pism
