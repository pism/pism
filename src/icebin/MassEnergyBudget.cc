#include <iostream>
#include "pism/icebin/MassEnergyBudget.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism{
namespace icebin{

MassEnthVec2S::MassEnthVec2S(std::shared_ptr<const Grid> my_grid,
                             const std::string &my_name)
  : mass(my_grid, my_name + ".mass"),
    enth(my_grid, my_name + ".enth")
{
  // empty
}

void MassEnthVec2S::set_attrs(const std::string &pism_intent, const std::string &long_name,
                              const std::string &units, const std::string &standard_name) {
  (void) standard_name;         // FIXME: these quantities don't have standard names
  auto mass_units   = "kg " + units;
  auto energy_units = "J " + units;
  mass.metadata().intent(pism_intent).long_name(long_name + " (mass portion)").units(mass_units);
  enth.metadata().intent(pism_intent).long_name(long_name + " (enthalpy portion)").units(energy_units);
}

MassEnergyBudget::MassEnergyBudget(std::shared_ptr<const Grid> grid, std::string const &prefix)
    : total(grid, prefix + "total"),
      basal_frictional_heating(grid, prefix + "basal_frictional_heating"),
      strain_heating(grid, prefix + "strain_heating"),
      geothermal_flux(grid, prefix + "geothermal_flux"),
      upward_geothermal_flux(grid, prefix + "upward_geothermal_flux"),
      calving(grid, prefix + "calving"),
      icebin_xfer(grid, prefix + "icebin_xfer"),
      icebin_deltah(grid, prefix + "icebin_deltah"),
      pism_smb(grid, prefix + "pism_smb"),
      href_to_h(grid, prefix + "href_to_h"),
      nonneg_rule(grid, prefix + "nonneg_rule"),
      melt_grounded(grid, prefix + "melt_grounded"),
      melt_floating(grid, prefix + "melt_floating"),
      internal_advection(grid, prefix + "internal_advection"),
      epsilon(grid, prefix + "epsilon") {
  printf("MassEnergyBudget(%p)::create()\n", (void *)this);

  // ----------- Mass and Enthalpy State of the Ice Sheet
  total.set_attrs("diagnostic", "State of the ice sheet (NOT a difference between timetseps)",
                  "m-2", "total");
  add_massenth(total, TOTAL, "", "");

  // ----------- Heat generation of flows [vertical]
  // Postive means heat is flowing INTO the ice sheet.
  basal_frictional_heating.metadata(0)
      .intent("internal")
      .long_name("Basal frictional heating")
      .units("W m-2")
      .glaciological_units("W m-2")
      .standard_name("");
  add_enth(basal_frictional_heating, DELTA, "basal_frictional_heating");

  strain_heating.metadata(0)
      .intent("internal")
      .long_name("Strain heating")
      .units("W m-2")
      .glaciological_units("W m-2")
      .standard_name("");
  add_enth(strain_heating, DELTA, "strain_heating"); //!< Total amount of strain heating [J/m^2]

  geothermal_flux.metadata(0)
      .intent("internal")
      .long_name("Geothermal energy through (compare to upward_geothermal_flux?)")
      .units("W m-2")
      .glaciological_units("W m-2")
      .standard_name("");
  add_enth(geothermal_flux, 0, "geothermal_flux"); //!< Total amount of geothermal energy [J/m^2]

  upward_geothermal_flux.metadata(0)
      .intent("internal")
      .long_name("Geothermal energy through (compare to geothermal_flux?)")
      .units("W m-2")
      .glaciological_units("W m-2")
      .standard_name("");
  add_enth(upward_geothermal_flux, DELTA,
           "upward_geothermal_flux"); //!< Total amount of geothermal energy [J/m^2]

  // ----------- Mass advection, with accompanying enthalpy change
  // Postive means mass/enthalpy is flowing INTO the ice sheet.
  std::string name;

  calving.set_attrs("diagnostic", "Mass/Enthalpy gain from calving.  Should be negative.",
                    "m-2 s-1", "calving");
  add_massenth(calving, DELTA, "calving.mass", "calving.enth");

  pism_smb.set_attrs("diagnostic", "pism_smb", "m-2 s-1", "pism_smb");
  // No DELTA< does not participate in epsilon computation
  add_massenth(pism_smb, 0, "pism_smb.mass", "pism_smb.enth");

  icebin_xfer.set_attrs("diagnostic", "icebin_xfer", "m-2 s-1", "icebin_xfer");
  add_massenth(icebin_xfer, DELTA, "icebin_xfer.mass", "icebin_xfer.enth");

  icebin_deltah.metadata(0)
      .intent("diagnostic")
      .long_name("icebin_deltah")
      .units("J m-2 s-1");
  add_enth(icebin_deltah, DELTA, "");

  href_to_h.metadata(0)
      .intent("diagnostic")
      .long_name("href_to_h")
      .units("kg m-2 s-1");
  add_mass(href_to_h, 0, "");

  nonneg_rule.metadata(0)
      .intent("diagnostic")
      .long_name("nonneg_rule")
      .units("kg m-2 s-1");
  add_mass(nonneg_rule, 0, "");


  melt_grounded.set_attrs("diagnostic", "Basal melting of grounded ice (negative)", "m-2 s-1",
                          "melt_grounded");
  add_massenth(melt_grounded, DELTA, "melt_grounded.mass", "melt_grounded.enth");

  melt_floating.set_attrs("diagnostic", "Sub-shelf melting (negative)", "m-2 s-1", "melt_floating");
  add_massenth(melt_floating, DELTA, "melt_floating.mass", "melt_floating.enth");

  // ----------- Advection WITHIN the ice sheet
  internal_advection.set_attrs("diagnostic", "Advection within the ice sheet", "m-2 s-1",
                               "internal_advection");
  add_massenth(internal_advection, DELTA, "internal_advection.mass", "internal_advection.enth");

  // ----------- Balance the Budget
  epsilon.set_attrs("diagnostic", "Unaccounted-for changes", "m-2 s-1", "epsilon");
  add_massenth(epsilon, EPSILON, "epsilon.mass", "epsilon.enth");
}

std::ostream &MassEnergyBudget::print_formulas(std::ostream &out)
{
	// MASS
	out << "epsilon.mass = total.mass -" << std::endl;
	out << "    (";
	for (auto ii=all_vecs.begin(); ii != all_vecs.end(); ++ii) {
		if ((ii->flags & (DELTA | MASS)) != (DELTA | MASS)) continue;
		char str[20];
		sprintf(str, "%p", (void*)&ii->vec);
		out << ii->vec.get_name() << " + ";
	}
	out << ")" << std::endl;

	// Energy
	out << "epsilon.enth = total.enth -" << std::endl;
	out << "    (";
	for (auto ii=all_vecs.begin(); ii != all_vecs.end(); ++ii) {
		if ((ii->flags & (DELTA | ENTH)) != (DELTA | ENTH)) continue;
		char str[20];
		sprintf(str, "%p", (void*)&ii->vec);
		out << ii->vec.get_name() << " + ";
	}
	out << ")" << std::endl;

	return out;
}


void MassEnergyBudget::set_epsilon(std::shared_ptr<const Grid> grid)
{
	// ==> epsilon = (sum of fluxes) - total
	printf("BEGIN MassEnergyBudget::set_epsilon()\n");

	// -------- Mass
	epsilon.mass.begin_access();
	total.mass.begin_access();
	for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
	for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
		epsilon.mass(i,j) = total.mass(i,j);
	}}
	total.mass.end_access();

	for (auto &ii : all_vecs) {
		if ((ii.flags & (DELTA | MASS)) != (DELTA | MASS)) continue;

		printf("epsilon.mass: %s\n", ii.vec.get_name().c_str());

		ii.vec.begin_access();
		for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
		for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
			epsilon.mass(i,j) -= ii.vec(i,j);
		}}
		ii.vec.end_access();
	}
	epsilon.mass.end_access();

	// -------- Energy
	epsilon.enth.begin_access();
	total.enth.begin_access();
	for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
	for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
		epsilon.enth(i,j) = total.enth(i,j);
	}}
	total.enth.end_access();

#if 1
	for (auto &ii : all_vecs) {
		if ((ii.flags & (DELTA | ENTH)) != (DELTA | ENTH)) continue;

		printf("epsilon.energy: %s\n", ii.vec.get_name().c_str());

		ii.vec.begin_access();
		for (int i = grid->xs(); i < grid->xs() + grid->xm(); ++i) {
		for (int j = grid->ys(); j < grid->ys() + grid->ym(); ++j) {
			epsilon.enth(i,j) -= ii.vec(i,j);
		}}
		ii.vec.end_access();
	}
#endif
	epsilon.enth.end_access();

	printf("END MassEnergyBudget::set_epsilon()\n");
}

}}
