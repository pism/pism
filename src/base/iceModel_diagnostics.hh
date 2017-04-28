// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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

namespace diagnostics {
//! \brief Computes vertically-averaged ice hardness.
class HardnessAverage : public Diag<IceModel>
{
public:
  HardnessAverage(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

class IceAreaFraction : public Diag<IceModel>
{
public:
  IceAreaFraction(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

class IceAreaFractionGrounded : public Diag<IceModel>
{
public:
  IceAreaFractionGrounded(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

class IceAreaFractionFloating : public Diag<IceModel>
{
public:
  IceAreaFractionFloating(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes a diagnostic field filled with processor rank values.
class Rank : public Diag<IceModel>
{
public:
  Rank(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes CTS, CTS = E/E_s(p).
class CTS : public Diag<IceModel>
{
public:
  CTS(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes ice temperature from enthalpy.
class Temperature : public Diag<IceModel>
{
public:
  Temperature(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature.
class TemperaturePA : public Diag<IceModel>
{
public:
  TemperaturePA(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes basal values of the pressure-adjusted temperature.
class TemperaturePABasal : public Diag<IceModel>
{
public:
  TemperaturePABasal(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes surface values of ice enthalpy.
class EnthalpySurface : public Diag<IceModel>
{
public:
  EnthalpySurface(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes enthalpy at the base of the ice.
class EnthalpyBasal : public Diag<IceModel>
{
public:
  EnthalpyBasal(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes ice temperature at the base of the ice.
class TemperatureBasal : public Diag<IceModel>
{
public:
  TemperatureBasal(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes ice temperature at the surface of the ice.
class TemperatureSurface : public Diag<IceModel>
{
public:
  TemperatureSurface(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes the liquid water fraction.
class LiquidFraction : public Diag<IceModel>
{
public:
  LiquidFraction(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes the total thickness of temperate ice in a column.
class TemperateIceThickness : public Diag<IceModel>
{
public:
  TemperateIceThickness(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};
//! \brief Computes the thickness of the basal layer of temperate ice.
class TemperateIceThicknessBasal : public Diag<IceModel>
{
public:
  TemperateIceThicknessBasal(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes the total ice volume in glacierized areas.
class VolumeGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total ice volume.
class VolumeNonGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeNonGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total ice volume, which is relevant for sea-level
class SeaLevelVolume : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  SeaLevelVolume(const IceModel *m);
  double compute();
};

//! \brief Computes the rate of change of the total ice volume in glacierized areas.
class VolumeRateOfChangeGlacierized : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  VolumeRateOfChangeGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the rate of change of the total ice volume.
class VolumeRateOfChangeNonGlacierized : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  VolumeRateOfChangeNonGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total ice area.
class AreaGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  AreaGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total ice mass in glacierized areas.
class MassGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  MassGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total ice mass.
class MassNonGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  MassNonGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total mass of the ice not displacing sea water.
class MassNotDisplacingSeaWater : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  MassNotDisplacingSeaWater(const IceModel *m);
  double compute();
};

//! \brief Computes the rate of change of the total ice mass in glacierized areas.
class MassRateOfChangeGlacierized : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  MassRateOfChangeGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the rate of change of the total ice mass.
class MassRateOfChangeNonGlacierized : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  MassRateOfChangeNonGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total volume of the temperate ice in glacierized areas.
class VolumeGlacierizedTemperate : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeGlacierizedTemperate(IceModel *m);
  double compute();
};

//! \brief Computes the total volume of the temperate ice.
class VolumeNonGlacierizedTemperate : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeNonGlacierizedTemperate(IceModel *m);
  double compute();
};

//! \brief Computes the total volume of the cold ice in glacierized areas.
class VolumeGlacierizedCold : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeGlacierizedCold(IceModel *m);
  double compute();
};

//! \brief Computes the total volume of the cold ice.
class VolumeNonGlacierizedCold : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeNonGlacierizedCold(IceModel *m);
  double compute();
};

//! \brief Computes the total area of the temperate ice.
class AreaGlacierizedTemperateBase : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  AreaGlacierizedTemperateBase(IceModel *m);
  double compute();
};

//! \brief Computes the total area of the cold ice.
class AreaGlacierizedColdBase : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  AreaGlacierizedColdBase(IceModel *m);
  double compute();
};

//! \brief Computes the total ice enthalpy in glacierized areas.
class EnthalpyGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  EnthalpyGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total ice enthalpy.
class EnthalpyNonGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  EnthalpyNonGlacierized(IceModel *m);
  double compute();
};

//! \brief Computes the total grounded ice area.
class AreaGlacierizedGrounded : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  AreaGlacierizedGrounded(IceModel *m);
  double compute();
};

//! \brief Computes the total floating ice area.
class AreaGlacierizedShelf : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  AreaGlacierizedShelf(IceModel *m);
  double compute();
};

//! \brief Computes the total grounded ice volume.
class VolumeGlacierizedGrounded : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeGlacierizedGrounded(IceModel *m);
  double compute();
};

//! \brief Computes the total floating ice volume.
class VolumeGlacierizedShelf : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  VolumeGlacierizedShelf(IceModel *m);
  double compute();
};

//! \brief Reports the mass continuity time step.
class TimeStepLength : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  TimeStepLength(const IceModel *m);
  double compute();
};

//! \brief Reports maximum diffusivity.
class MaxDiffusivity : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  MaxDiffusivity(const IceModel *m);
  double compute();
};

//! \brief Reports the total surface ice flux.
class MassFluxSurface : public TSDiag<TSFluxDiagnostic, IceModel> // FIXME_
{
public:
  MassFluxSurface(const IceModel *m);
  double compute();
};

//! \brief Reports the total basal ice flux over the grounded region.
class MassFluxBasalGrounded : public TSDiag<TSFluxDiagnostic, IceModel> // FIXME_
{
public:
  MassFluxBasalGrounded(const IceModel *m);
  double compute();
};

//! \brief Reports the total sub-shelf ice flux.
class MassFluxBasalFloating : public TSDiag<TSFluxDiagnostic, IceModel> // FIXME_
{
public:
  MassFluxBasalFloating(const IceModel *m);
  double compute();
};

//! \brief Reports the total discharge flux.
class MassFluxDischarge : public TSDiag<TSFluxDiagnostic, IceModel> // FIXME_
{
public:
  MassFluxDischarge(const IceModel *m);
  double compute();
};

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
class MaxHorizontalVelocity : public TSDiag<TSSnapshotDiagnostic, IceModel> // GOOD
{
public:
  MaxHorizontalVelocity(const IceModel *m);
  double compute();
};

//! \brief Computes the 2D height above flotation.
class HeightAboveFloatation : public Diag<IceModel>
{
public:
  HeightAboveFloatation(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};


//! \brief Computes the mass per cell.
class IceMass : public Diag<IceModel>
{
public:
  IceMass(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

/*! @brief Sea-level adjusted bed topography (zero at sea level). */
class BedTopographySeaLevelAdjusted : public Diag<IceModel>
{
public:
  BedTopographySeaLevelAdjusted(const IceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Ice hardness computed using the SIA flow law. */
class Hardness : public Diag<IceModel>
{
public:
  Hardness(const IceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Effective viscosity of ice (3D). */
class Viscosity : public Diag<IceModel>
{
public:
  Viscosity(IceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes latitude and longitude bounds.
class LatLonBounds : public Diag<IceModel>
{
public:
  LatLonBounds(const IceModel *m,
                          const std::string &var_name,
                          const std::string &proj_string);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
protected:
  std::string m_var_name, m_proj_string;
};

} // end of namespace diagnostics
} // end of namespace pism

#endif  /* _ICEMODEL_DIAGNOSTICS_H_ */
