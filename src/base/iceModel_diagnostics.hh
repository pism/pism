// Copyright (C) 2010, 2011, 2012 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef _ICEMODEL_DIAGNOSTICS_H_
#define _ICEMODEL_DIAGNOSTICS_H_

#include "iceModel.hh"
#include "PISMDiagnostic.hh"

//! \brief Computes vertically-averaged ice hardness.
class IceModel_hardav : public PISMDiag<IceModel>
{
public:
  IceModel_hardav(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes a diagnostic field filled with processor rank values.
class IceModel_rank : public PISMDiag<IceModel>
{
public:
  IceModel_rank(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes CTS, CTS = E/E_s(p).
class IceModel_cts : public PISMDiag<IceModel>
{
public:
  IceModel_cts(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the number of ice-filled cells is a processor's domain.
class IceModel_proc_ice_area : public PISMDiag<IceModel>
{
public:
  IceModel_proc_ice_area(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ice temperature from enthalpy.
class IceModel_temp : public PISMDiag<IceModel>
{
public:
  IceModel_temp(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature.
class IceModel_temp_pa : public PISMDiag<IceModel>
{
public:
  IceModel_temp_pa(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes basal values of the pressure-adjusted temperature.
class IceModel_temppabase : public PISMDiag<IceModel>
{
public:
  IceModel_temppabase(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes surface values of ice enthalpy.
class IceModel_enthalpysurf : public PISMDiag<IceModel>
{
public:
  IceModel_enthalpysurf(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes enthalpy at the base of the ice.
class IceModel_enthalpybase : public PISMDiag<IceModel>
{
public:
  IceModel_enthalpybase(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ice temperature at the base of the ice.
class IceModel_tempbase : public PISMDiag<IceModel>
{
public:
  IceModel_tempbase(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ice temperature at the surface of the ice.
class IceModel_tempsurf : public PISMDiag<IceModel>
{
public:
  IceModel_tempsurf(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the liquid water fraction.
class IceModel_liqfrac : public PISMDiag<IceModel>
{
public:
  IceModel_liqfrac(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the total thickness of temperate ice in a column.
class IceModel_tempicethk : public PISMDiag<IceModel>
{
public:
  IceModel_tempicethk(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};
//! \brief Computes the thickness of the basal layer of temperate ice.
class IceModel_tempicethk_basal : public PISMDiag<IceModel>
{
public:
  IceModel_tempicethk_basal(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the total ice volume.
class IceModel_ivol : public PISMTSDiag<IceModel>
{
public:
  IceModel_ivol(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total ice volume, which is relevant for sea-level
class IceModel_slvol : public PISMTSDiag<IceModel>
{
public:
  IceModel_slvol(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the rate of change of the total ice volume.
class IceModel_divoldt : public PISMTSDiag<IceModel>
{
public:
  IceModel_divoldt(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total ice area.
class IceModel_iarea : public PISMTSDiag<IceModel>
{
public:
  IceModel_iarea(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total ice mass.
class IceModel_imass : public PISMTSDiag<IceModel>
{
public:
  IceModel_imass(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the rate of change of the total ice mass.
class IceModel_dimassdt : public PISMTSDiag<IceModel>
{
public:
  IceModel_dimassdt(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total volume of the temperate ice.
class IceModel_ivoltemp : public PISMTSDiag<IceModel>
{
public:
  IceModel_ivoltemp(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total volume of the cold ice.
class IceModel_ivolcold : public PISMTSDiag<IceModel>
{
public:
  IceModel_ivolcold(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total area of the temperate ice.
class IceModel_iareatemp : public PISMTSDiag<IceModel>
{
public:
  IceModel_iareatemp(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total area of the cold ice.
class IceModel_iareacold : public PISMTSDiag<IceModel>
{
public:
  IceModel_iareacold(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total ice enthalpy.
class IceModel_ienthalpy : public PISMTSDiag<IceModel>
{
public:
  IceModel_ienthalpy(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total grounded ice area.
class IceModel_iareag : public PISMTSDiag<IceModel>
{
public:
  IceModel_iareag(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total floating ice area.
class IceModel_iareaf : public PISMTSDiag<IceModel>
{
public:
  IceModel_iareaf(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total grounded ice volume.
class IceModel_ivolg : public PISMTSDiag<IceModel>
{
public:
  IceModel_ivolg(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Computes the total floating ice volume.
class IceModel_ivolf : public PISMTSDiag<IceModel>
{
public:
  IceModel_ivolf(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the mass continuity time step.
class IceModel_dt : public PISMTSDiag<IceModel>
{
public:
  IceModel_dt(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports maximum diffusivity.
class IceModel_max_diffusivity : public PISMTSDiag<IceModel>
{
public:
  IceModel_max_diffusivity(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the total surface ice flux.
class IceModel_surface_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_surface_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the cumulative total surface ice flux.
class IceModel_surface_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_surface_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the total basal ice flux over the grounded region.
class IceModel_grounded_basal_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the cumulative total basal ice flux over the grounded region.
class IceModel_grounded_basal_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the total sub-shelf ice flux.
class IceModel_sub_shelf_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_sub_shelf_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the total sub-shelf ice flux.
class IceModel_sub_shelf_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_sub_shelf_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the 'numerical' ice flux resulting from enforcing the 'thk
//! >= 0' rule.
class IceModel_nonneg_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_nonneg_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the cumulative 'numerical' ice flux resulting from enforcing the 'thk
//! >= 0' rule.
class IceModel_nonneg_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_nonneg_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the -ocean_kill flux.
class IceModel_ocean_kill_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_ocean_kill_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the cumulative -ocean_kill flux.
class IceModel_ocean_kill_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_ocean_kill_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the total -float_kill flux.
class IceModel_float_kill_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_float_kill_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the cumulative total -float_kill flux.
class IceModel_float_kill_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_float_kill_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the total discharge flux.
class IceModel_discharge_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_discharge_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the cumulative total discharge flux.
class IceModel_discharge_flux_cumulative : public PISMTSDiag<IceModel>
{
public:
  IceModel_discharge_flux_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports cumulative surface mass balance.
class IceModel_climatic_mass_balance_cumulative : public PISMDiag<IceModel>
{
public:
  IceModel_climatic_mass_balance_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Reports cumulative ocean kill flux.
class IceModel_ocean_kill_flux_2D_cumulative : public PISMDiag<IceModel>
{
public:
  IceModel_ocean_kill_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes ocean_kill_flux, the calving flux due to the "-ocean_kill"
//! mechanism.
class IceModel_ocean_kill_flux_2D : public PISMDiag<IceModel>
{
public:
  IceModel_ocean_kill_flux_2D(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
  virtual PetscErrorCode update_cumulative();
protected:
  IceModelVec2S last_ocean_kill_flux_cumulative;
  PetscReal last_report_time;
};

//! \brief Computes dHdt, the ice thickness rate of change.
class IceModel_dHdt : public PISMDiag<IceModel>
{
public:
  IceModel_dHdt(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
  virtual PetscErrorCode update_cumulative();
protected:
  IceModelVec2S last_ice_thickness;
  PetscReal last_report_time;
};

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
class IceModel_max_hor_vel : public PISMTSDiag<IceModel>
{
public:
  IceModel_max_hor_vel(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the mass flux from the mass tracked using ice thickness
//! (thk) to the mass tracked using Href.
class IceModel_H_to_Href_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_H_to_Href_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the mass flux from the mass tracked using Href to the mass
//! tracked using ice thickness (thk).
class IceModel_Href_to_H_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_Href_to_H_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the sum(div Q) flux (to diagnose issues in the mass
//! transport scheme).
class IceModel_sum_divQ_flux : public PISMTSDiag<IceModel>
{
public:
  IceModel_sum_divQ_flux(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode update(PetscReal a, PetscReal b);
};

//! \brief Reports the 2D cumulative (numerical) flux due to enforcing
//! non-negativity of ice thickness.
class IceModel_nonneg_flux_2D_cumulative : public PISMDiag<IceModel>
{
public:
  IceModel_nonneg_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


//! \brief Reports the 2D cumulative grounded basal flux.
class IceModel_grounded_basal_flux_2D_cumulative : public PISMDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Reports the 2D cumulative floating basal flux.
class IceModel_floating_basal_flux_2D_cumulative : public PISMDiag<IceModel>
{
public:
  IceModel_floating_basal_flux_2D_cumulative(IceModel *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};


#endif  /* _ICEMODEL_DIAGNOSTICS_H_ */

