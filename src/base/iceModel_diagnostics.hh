// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include "base/util/PISMDiagnostic.hh"

namespace pism {

//! \brief Computes vertically-averaged ice hardness.
class IceModel_hardav : public Diag<IceModel>
{
public:
  IceModel_hardav(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

class IceModel_land_ice_area_fraction : public Diag<IceModel>
{
public:
  IceModel_land_ice_area_fraction(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

class IceModel_grounded_ice_sheet_area_fraction : public Diag<IceModel>
{
public:
  IceModel_grounded_ice_sheet_area_fraction(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

class IceModel_floating_ice_sheet_area_fraction : public Diag<IceModel>
{
public:
  IceModel_floating_ice_sheet_area_fraction(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes a diagnostic field filled with processor rank values.
class IceModel_rank : public Diag<IceModel>
{
public:
  IceModel_rank(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes CTS, CTS = E/E_s(p).
class IceModel_cts : public Diag<IceModel>
{
public:
  IceModel_cts(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes the number of ice-filled cells is a processor's domain.
class IceModel_proc_ice_area : public Diag<IceModel>
{
public:
  IceModel_proc_ice_area(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes ice temperature from enthalpy.
class IceModel_temp : public Diag<IceModel>
{
public:
  IceModel_temp(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature.
class IceModel_temp_pa : public Diag<IceModel>
{
public:
  IceModel_temp_pa(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes basal values of the pressure-adjusted temperature.
class IceModel_temppabase : public Diag<IceModel>
{
public:
  IceModel_temppabase(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes surface values of ice enthalpy.
class IceModel_enthalpysurf : public Diag<IceModel>
{
public:
  IceModel_enthalpysurf(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes enthalpy at the base of the ice.
class IceModel_enthalpybase : public Diag<IceModel>
{
public:
  IceModel_enthalpybase(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes ice temperature at the base of the ice.
class IceModel_tempbase : public Diag<IceModel>
{
public:
  IceModel_tempbase(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes ice temperature at the surface of the ice.
class IceModel_tempsurf : public Diag<IceModel>
{
public:
  IceModel_tempsurf(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes the liquid water fraction.
class IceModel_liqfrac : public Diag<IceModel>
{
public:
  IceModel_liqfrac(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes the total thickness of temperate ice in a column.
class IceModel_tempicethk : public Diag<IceModel>
{
public:
  IceModel_tempicethk(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};
//! \brief Computes the thickness of the basal layer of temperate ice.
class IceModel_tempicethk_basal : public Diag<IceModel>
{
public:
  IceModel_tempicethk_basal(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};
//! \brief Computes the flux divergence.
class IceModel_flux_divergence : public Diag<IceModel>
{
public:
  IceModel_flux_divergence(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes the total ice volume in glacierized areas.
class IceModel_volume_glacierized : public TSDiag<IceModel>
{
public:
  IceModel_volume_glacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice volume.
class IceModel_volume_nonglacierized : public TSDiag<IceModel>
{
public:
  IceModel_volume_nonglacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice volume, which is relevant for sea-level
class IceModel_slvol : public TSDiag<IceModel>
{
public:
  IceModel_slvol(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the rate of change of the total ice volume in glacierized areas.
class IceModel_volume_rate_of_change_glacierized : public TSDiag<IceModel>
{
public:
  IceModel_volume_rate_of_change_glacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the rate of change of the total ice volume.
class IceModel_volume_rate_of_change_nonglacierized : public TSDiag<IceModel>
{
public:
  IceModel_volume_rate_of_change_nonglacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice area.
class IceModel_area_glacierized : public TSDiag<IceModel>
{
public:
  IceModel_area_glacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice mass in glacierized areas.
class IceModel_mass_glacierized : public TSDiag<IceModel>
{
public:
  IceModel_mass_glacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice mass.
class IceModel_mass_nonglacierized : public TSDiag<IceModel>
{
public:
  IceModel_mass_nonglacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total mass of the ice not displacing sea water.
class IceModel_limnsw : public TSDiag<IceModel>
{
public:
  IceModel_limnsw(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the rate of change of the total ice mass in glacierized areas.
class IceModel_mass_rate_of_change_glacierized : public TSDiag<IceModel>
{
public:
  IceModel_mass_rate_of_change_glacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the rate of change of the total ice mass.
class IceModel_mass_rate_of_change_nonglacierized : public TSDiag<IceModel>
{
public:
  IceModel_mass_rate_of_change_nonglacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total volume of the temperate ice in glacierized areas.
class IceModel_volume_glacierized_temperate : public TSDiag<IceModel>
{
public:
  IceModel_volume_glacierized_temperate(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total volume of the temperate ice.
class IceModel_volume_nonglacierized_temperate : public TSDiag<IceModel>
{
public:
  IceModel_volume_nonglacierized_temperate(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total volume of the cold ice in glacierized areas.
class IceModel_volume_glacierized_cold : public TSDiag<IceModel>
{
public:
  IceModel_volume_glacierized_cold(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total volume of the cold ice.
class IceModel_volume_nonglacierized_cold : public TSDiag<IceModel>
{
public:
  IceModel_volume_nonglacierized_cold(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total area of the temperate ice.
class IceModel_area_glacierized_temperate_base : public TSDiag<IceModel>
{
public:
  IceModel_area_glacierized_temperate_base(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total area of the cold ice.
class IceModel_area_glacierized_cold_base : public TSDiag<IceModel>
{
public:
  IceModel_area_glacierized_cold_base(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice enthalpy in glacierized areas.
class IceModel_enthalpy_glacierized : public TSDiag<IceModel>
{
public:
  IceModel_enthalpy_glacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice enthalpy.
class IceModel_enthalpy_nonglacierized : public TSDiag<IceModel>
{
public:
  IceModel_enthalpy_nonglacierized(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total grounded ice area.
class IceModel_area_glacierized_grounded : public TSDiag<IceModel>
{
public:
  IceModel_area_glacierized_grounded(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total floating ice area.
class IceModel_area_glacierized_shelf : public TSDiag<IceModel>
{
public:
  IceModel_area_glacierized_shelf(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total grounded ice volume.
class IceModel_volume_glacierized_grounded : public TSDiag<IceModel>
{
public:
  IceModel_volume_glacierized_grounded(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total floating ice volume.
class IceModel_volume_glacierized_shelf : public TSDiag<IceModel>
{
public:
  IceModel_volume_glacierized_shelf(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the mass continuity time step.
class IceModel_dt : public TSDiag<IceModel>
{
public:
  IceModel_dt(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports maximum diffusivity.
class IceModel_max_diffusivity : public TSDiag<IceModel>
{
public:
  IceModel_max_diffusivity(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total surface ice flux.
class IceModel_surface_flux : public TSDiag<IceModel>
{
public:
  IceModel_surface_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative total surface ice flux.
class IceModel_surface_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_surface_flux_cumulative(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total basal ice flux over the grounded region.
class IceModel_grounded_basal_flux : public TSDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative total basal ice flux over the grounded region.
class IceModel_grounded_basal_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux_cumulative(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total sub-shelf ice flux.
class IceModel_sub_shelf_flux : public TSDiag<IceModel>
{
public:
  IceModel_sub_shelf_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total sub-shelf ice flux.
class IceModel_sub_shelf_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_sub_shelf_flux_cumulative(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the 'numerical' ice flux resulting from enforcing the 'thk
//! >= 0' rule.
class IceModel_nonneg_flux : public TSDiag<IceModel>
{
public:
  IceModel_nonneg_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative 'numerical' ice flux resulting from enforcing the 'thk
//! >= 0' rule.
class IceModel_nonneg_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_nonneg_flux_cumulative(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total discharge flux.
class IceModel_discharge_flux : public TSDiag<IceModel>
{
public:
  IceModel_discharge_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative total discharge flux.
class IceModel_discharge_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_discharge_flux_cumulative(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports cumulative surface mass balance.
class IceModel_climatic_mass_balance_cumulative : public Diag<IceModel>
{
public:
  IceModel_climatic_mass_balance_cumulative(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Computes dHdt, the ice thickness rate of change.
class IceModel_dHdt : public Diag<IceModel>
{
public:
  IceModel_dHdt(const IceModel *m);
  virtual void update_cumulative();
protected:
  virtual IceModelVec::Ptr compute_impl();
protected:
  IceModelVec2S m_last_ice_thickness;
  double m_last_report_time;
};

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
class IceModel_max_hor_vel : public TSDiag<IceModel>
{
public:
  IceModel_max_hor_vel(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the mass flux from the mass tracked using ice thickness
//! (thk) to the mass tracked using Href.
class IceModel_H_to_Href_flux : public TSDiag<IceModel>
{
public:
  IceModel_H_to_Href_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the mass flux from the mass tracked using Href to the mass
//! tracked using ice thickness (thk).
class IceModel_Href_to_H_flux : public TSDiag<IceModel>
{
public:
  IceModel_Href_to_H_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the sum(div Q) flux (to diagnose issues in the mass
//! transport scheme).
class IceModel_sum_divQ_flux : public TSDiag<IceModel>
{
public:
  IceModel_sum_divQ_flux(const IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the 2D cumulative (numerical) flux due to enforcing
//! non-negativity of ice thickness.
class IceModel_nonneg_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_nonneg_flux_2D_cumulative(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Reports the 2D cumulative grounded basal flux.
class IceModel_grounded_basal_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_grounded_basal_flux_2D_cumulative(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Reports the 2D cumulative floating basal flux.
class IceModel_floating_basal_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_floating_basal_flux_2D_cumulative(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Reports the 2D cumulative discharge (calving) flux.
class IceModel_discharge_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_discharge_flux_2D_cumulative(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

//! \brief Reports the 2D cumulative discharge (calving) flux.
class IceModel_discharge_flux_2D : public Diag<IceModel>
{
public:
  IceModel_discharge_flux_2D(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
  IceModelVec2S m_last_cumulative_discharge;
  double m_last_report_time;
};

//! \brief Reports surface mass balance flux (average over reporting interval).
class IceModel_surface_mass_balance_average : public Diag<IceModel>
{
public:
  IceModel_surface_mass_balance_average(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
  IceModelVec2S m_last_cumulative_SMB;
  double m_last_report_time;
};

//! \brief Reports the 2D cumulative discharge (calving) flux.
class IceModel_basal_mass_balance_average : public Diag<IceModel>
{
public:
  IceModel_basal_mass_balance_average(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
  IceModelVec2S m_last_cumulative_BMB;
  double m_last_report_time;
};

//! \brief Computes the 2D height above flotation.
class IceModel_height_above_flotation : public Diag<IceModel>
{
public:
  IceModel_height_above_flotation(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};


//! \brief Computes the mass per cell.
class IceModel_ice_mass : public Diag<IceModel>
{
public:
  IceModel_ice_mass(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl();
};

/*! @brief Sea-level adjusted bed topography (zero at sea level). */
class IceModel_topg_sl_adjusted : public Diag<IceModel>
{
public:
  IceModel_topg_sl_adjusted(const IceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Ice hardness computed using the SIA flow law. */
class IceModel_hardness : public Diag<IceModel>
{
public:
  IceModel_hardness(const IceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Effective viscosity of ice (3D). */
class IceModel_viscosity : public Diag<IceModel>
{
public:
  IceModel_viscosity(IceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

} // end of namespace pism

#if (PISM_USE_PROJ4==1)
#include <proj_api.h>
#endif

namespace pism {

//! \brief Computes latitude and longitude bounds.
class IceModel_lat_lon_bounds : public Diag<IceModel>
{
public:
  IceModel_lat_lon_bounds(const IceModel *m,
                          const std::string &var_name,
                          const std::string &proj_string);
protected:
  virtual IceModelVec::Ptr compute_impl();
protected:
  std::string m_var_name, m_proj_string;
};

} // end of namespace pism

#endif  /* _ICEMODEL_DIAGNOSTICS_H_ */
