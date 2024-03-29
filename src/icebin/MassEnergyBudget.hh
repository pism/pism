#pragma once

#include "pism/util/array/Scalar.hh"
#include <memory>

namespace pism {
namespace icebin {

/** Encapsulates mass and enthalpy together.  Used to tabulate total
enthalpy of a bunch of advected H2O based on its mass and specific
enthalpy.  This allows us to have only one C++ variable per advected
quantity, instead of two. */
struct MassEnthVec2S : public pism::PetscAccessible {
  pism::array::Scalar mass;
  pism::array::Scalar enth;

  MassEnthVec2S(std::shared_ptr<const pism::Grid> grid, const std::string &name);

  void set_attrs(const std::string &long_name, const std::string &units);

  virtual void begin_access() const {
    mass.begin_access();
    enth.begin_access();
  }

  virtual void end_access() const {
    mass.end_access();
    enth.end_access();
  }

#if 0
    /** @param mass [kg m-2]
    @param specific_enthalpy [J kg-1] */
    void add_mass(int i, int j, double mass, double specific_enthalpy) {
        mass(i,j) += mass;
        enthalpy(i,j) += mass * specific_enthalpy;
    }

    void clear() {
        mass.set(0);
        enth.set(0);
    }
#endif
};

struct VecWithFlags {
  pism::array::Scalar &vec;
  int flags;

  /** IF this variable is used directly as a contractual ice model
    output, the name of that contract entry. */
  std::string contract_name;

  VecWithFlags(pism::array::Scalar &_vec, int _flags, std::string const &_contract_name)
      : vec(_vec), flags(_flags), contract_name(_contract_name) {
  }
};

class MassEnergyBudget {
public:
  // ============================================================
  // Total State

  // ------------ Enthalpy State
  MassEnthVec2S total; // Total mass [kg m-2] and enthalpy [J m-2] of the ice sheet

  // =============================================================
  // Cumulative Fluxes

  // ======================= Variables to accumulate PISM output
  // These are accumulated as [kg m-2] or [J m-2]
  // They are accumulated for the life of the simulation, and never zeroed.
  // Other instances of MassEnergyBudget are used to compute differences.

  // ----------- Heat generation of flows [vertical]
  pism::array::Scalar basal_frictional_heating; //!< Total amount of basal friction heating [J/m^2]
  pism::array::Scalar strain_heating;           //!< Total amount of strain heating [J/m^2]
  pism::array::Scalar geothermal_flux;          //!< Total amount of geothermal energy [J/m^2]
  pism::array::Scalar upward_geothermal_flux;   //!< Total amount of geothermal energy [J/m^2]

  // ----------- Mass advection, with accompanying enthalpy change
  // The enthalpy reported for these mass fluxes are the enthalpies
  // AS REPORTED TO ICEBIN!  That is not necessarily the same as the enthalpy
  // that PISM sees internally.
  MassEnthVec2S calving; //!< Equal to IceModel::discharge_flux_2D_cumulative

  // Let's not use basal_runoff (which would be derived from the variables
  // in NullTransportHydrology).
  //
  // 1. It is derived directly from
  //    bmelt/basal_meltrate/melt_grounded/meltrate_grounded, which is
  //    already computed here.
  //
  // 2. bmelt plays directly into removal of mass from the ice sheet (see
  //    accumulateFluxes_massContExplicitStep() in iMgeometry.cc).
  //    Including basal_runoff would be double-counting it.
  //  MassEnthVec2S basal_runoff;     //!< Enthalpy here is predictable, since runoff is 0C 100% water fraction.

  MassEnthVec2S smb;          //!< accumulation / ablation, as provided by Icebin
  pism::array::Scalar deltah; //!< Change in enthalpy of top layer
  MassEnthVec2S
      pism_smb; //! SMB as seen by PISM in iMgeometry.cc massContExplicitSte().  Used to check icebin_smb.mass, but does not figure into contract.
  pism::array::Scalar href_to_h;
  pism::array::Scalar nonneg_rule;
  MassEnthVec2S melt_grounded; //!< basal melt (grounded) (from summing meltrate_grounded)
  MassEnthVec2S melt_floating; //!< sub-shelf melt (from summing meltrate_floating)

  // ----------- Mass advection WITHIN the ice sheet
  MassEnthVec2S internal_advection;
  //  MassEnthVec2S divQ_SIA;
  //  MassEnthVec2S divQ_SSA;     //!< flux divergence

  // ======================= Balance the Budget
  // At each step, we set epsilon as follows:
  // total - (sum of fluxes) + epsilon = 0
  // ==> epsilon = (sum of fluxes) - total
  MassEnthVec2S epsilon;

  // ======================== Different sets (flags)
  static const int MASS = 1;
  static const int ENTH = 2;

  // No post-processing: Value written to NetCDF is same as accumulated value.
  // Good for things taht are set once.
  // (eg: total.mass and total.enth, which are not wrten to NetCDF).
  static const int TOTAL = 4;
  // Variable to be accumulated, differenced and /dt at end;
  // and also go into computation of epsilon
  static const int DELTA = 8;
  // Do NOT include in computation of epsilon, since this flag IS
  // for the varaiable epsilon (and also pism_smb)
  static const int EPSILON = 16; // To be differenced at the end.
  static const int ADVECTION =
      32; // This energy term is due to advection of a related mass term (not acted upon)

  // ======================== Summary of above variables
  // This makes it easy to difference two MassEnergyBudget instances.
  std::vector<VecWithFlags> all_vecs;


  // =====================================================================
  std::ostream &print_formulas(std::ostream &out);

protected:
  void add_mass(pism::array::Scalar &vec, int flags, std::string const &contract_name) {
    all_vecs.push_back(VecWithFlags(vec, MASS | flags, contract_name));
  }

  /** @param contract_name The name of this variable in the ice model's output contract. */
  void add_enth(pism::array::Scalar &vec, int flags, std::string const &contract_name) {
    all_vecs.push_back(VecWithFlags(vec, ENTH | flags, contract_name));
  }

  void add_massenth(MassEnthVec2S &massenth, int flags, std::string const &contract_name_mass,
                    std::string const &contract_name_enth) {
    all_vecs.push_back(VecWithFlags(massenth.mass, MASS | flags, contract_name_mass));
    all_vecs.push_back(VecWithFlags(massenth.enth, ENTH | ADVECTION | flags, contract_name_enth));
  }

public:
  MassEnergyBudget(std::shared_ptr<const pism::Grid> grid, std::string const &prefix);

  void set_epsilon();
};

} // end of namespace icebin
} // end of namespace pism
