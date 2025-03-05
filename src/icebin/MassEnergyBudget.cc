#include <iostream>

#include "pism/icebin/MassEnergyBudget.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {
namespace icebin {

MassEnthVec2S::MassEnthVec2S(std::shared_ptr<const pism::Grid> grid, const std::string &name)
    : mass(grid, name + ".mass"), enth(grid, name + ".enth") {
}

void MassEnthVec2S::set_attrs(const std::string &long_name, const std::string &units) {
  mass.metadata().long_name(long_name + " (mass portion)").units("kg " + units);
  enth.metadata().long_name(long_name + " (enthalpy portion)").units("J " + units);
}

MassEnergyBudget::MassEnergyBudget(std::shared_ptr<const pism::Grid> grid,
                                   std::string const &prefix)
    : total(grid, prefix + "total"),
      basal_frictional_heating(grid, prefix + "basal_frictional_heating"),
      strain_heating(grid, prefix + "strain_heating"),
      geothermal_flux(grid, prefix + "geothermal_flux"),
      upward_geothermal_flux(grid, prefix + "upward_geothermal_flux"),
      calving(grid, prefix + "calving"),
      smb(grid, prefix + "smb"),
      deltah(grid, prefix + "deltah"),
      pism_smb(grid, prefix + "pism_smb"),
      href_to_h(grid, prefix + "href_to_h"),
      nonneg_rule(grid, prefix + "nonneg_rule"),
      melt_grounded(grid, prefix + "melt_grounded"),
      melt_floating(grid, prefix + "melt_floating"),
      internal_advection(grid, prefix + "internal_advection"),
      epsilon(grid, prefix + "epsilon") {

  // ----------- Mass and Enthalpy State of the Ice Sheet

  total.set_attrs("State of the ice sheet (NOT a difference between time steps)", "m-2");
  add_massenth(total, TOTAL, "", "");

  // ----------- Heat generation of flows [vertical]
  // Postive means heat is flowing INTO the ice sheet.

  basal_frictional_heating.metadata(0).long_name("Basal frictional heating").units("W m^-2");
  add_enth(basal_frictional_heating, DELTA, "basal_frictional_heating");

  strain_heating.metadata(0).long_name("Strain heating").units("W m^-2");
  add_enth(strain_heating, DELTA, "strain_heating"); //!< Total amount of strain heating [J/m^2]

  geothermal_flux.metadata(0)
      .long_name("Geothermal energy through (compare to upward_geothermal_flux?)")
      .units("W m^-2");
  add_enth(geothermal_flux, 0, "geothermal_flux"); //!< Total amount of geothermal energy [J/m^2]

  upward_geothermal_flux.metadata(0)
      .long_name("Geothermal energy through (compare to geothermal_flux?)")
      .units("W m^-2");
  // Total amount of geothermal energy [J/m^2]
  add_enth(upward_geothermal_flux, DELTA, "upward_geothermal_flux");

  // ----------- Mass advection, with accompanying enthalpy change
  // Postive means mass/enthalpy is flowing INTO the ice sheet.
  calving.set_attrs("Mass/Enthalpy gain from calving.  Should be negative.",
                    "m-2 s-1");
  add_massenth(calving, DELTA, "calving.mass", "calving.enth");

  // SMB as seen by PISM in iMgeometry.cc massContExplicitSte().
  // Used to check icebin_smb.mass, but does not figure into
  // contract.
  pism_smb.set_attrs("pism_smb", "m-2 s-1");
  // No DELTA< does not participate in epsilon computation
  add_massenth(pism_smb, EPSILON, "pism_smb.mass", "pism_smb.enth");

  // accumulation / ablation, as provided by Icebin
  smb.set_attrs("smb", "m-2 s-1");
  add_massenth(smb, DELTA, "smb.mass", "smb.enth");

  deltah.metadata(0).long_name("deltah").units("J m^-2 s^-1");
  add_enth(deltah, DELTA, "");

  href_to_h.metadata(0).long_name("href_to_h").units("kg m^-2 s^-1");
  add_mass(href_to_h, 0, "");

  nonneg_rule.metadata(0).long_name("nonneg_rule").units("kg m^-2 s^-1");
  add_mass(nonneg_rule, 0, "");


  melt_grounded.set_attrs("Basal melting of grounded ice (negative)", "m-2 s-1");
  add_massenth(melt_grounded, DELTA, "melt_grounded.mass", "melt_grounded.enth");

  melt_floating.set_attrs("Sub-shelf melting (negative)", "m-2 s-1");
  add_massenth(melt_floating, DELTA, "melt_floating.mass", "melt_floating.enth");

  // ----------- Advection WITHIN the ice sheet
  internal_advection.set_attrs("Advection within the ice sheet", "m-2 s-1");
  add_massenth(internal_advection, DELTA, "internal_advection.mass", "internal_advection.enth");

  // ----------- Balance the Budget
  epsilon.set_attrs("Unaccounted-for changes", "m-2 s-1");
  add_massenth(epsilon, EPSILON, "epsilon.mass", "epsilon.enth");
}

std::ostream &MassEnergyBudget::print_formulas(std::ostream &out) {
  // MASS
  out << "epsilon.mass = total.mass -" << std::endl;
  out << "    (";
  for (auto ii = all_vecs.begin(); ii != all_vecs.end(); ++ii) {
    if ((ii->flags & (DELTA | MASS)) != (DELTA | MASS)) {
      continue;
    }
    char str[20];
    sprintf(str, "%p", (void *)&ii->vec);
    out << ii->vec.get_name() << " + ";
  }
  out << ")" << std::endl;

  // Energy
  out << "epsilon.enth = total.enth -" << std::endl;
  out << "    (";
  for (auto ii = all_vecs.begin(); ii != all_vecs.end(); ++ii) {
    if ((ii->flags & (DELTA | ENTH)) != (DELTA | ENTH)) {
      continue;
    }
    char str[20];
    sprintf(str, "%p", (void *)&ii->vec);
    out << ii->vec.get_name() << " + ";
  }
  out << ")" << std::endl;

  return out;
}


void MassEnergyBudget::set_epsilon() {
  auto grid = epsilon.mass.grid();
  // ==> epsilon = (sum of fluxes) - total

  // -------- Mass
  epsilon.mass.begin_access();
  total.mass.begin_access();
  for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
    for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
      epsilon.mass(i, j) = total.mass(i, j);
    }
  }
  total.mass.end_access();

  for (auto &ii : all_vecs) {
    // This normally does not include things marked EPSILON
    // (which is mutually exclusive with DELTA)
    // This excluded epsilon (ourself), as well as redundant pism_smb
    if ((ii.flags & (DELTA | MASS)) != (DELTA | MASS)) {
      continue;
    }

    ii.vec.begin_access();
    for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
      for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
        epsilon.mass(i, j) -= ii.vec(i, j);
      }
    }
    ii.vec.end_access();
  }
  epsilon.mass.end_access();

  // -------- Energy
  epsilon.enth.begin_access();
  total.enth.begin_access();
  for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
    for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
      epsilon.enth(i, j) = total.enth(i, j);
    }
  }
  total.enth.end_access();

#if 1
  for (auto &ii : all_vecs) {
    if ((ii.flags & (DELTA | ENTH)) != (DELTA | ENTH)) {
      continue;
    }

    printf("epsilon.energy: %s\n", ii.vec.get_name().c_str());

    ii.vec.begin_access();
    for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
      for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
        epsilon.enth(i, j) -= ii.vec(i, j);
      }
    }
    ii.vec.end_access();
  }
#endif
  epsilon.enth.end_access();
}

} // namespace icebin
} // namespace pism
