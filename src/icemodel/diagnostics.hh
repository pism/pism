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

#include "IceModel.hh"
#include "pism/util/Diagnostic.hh"

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
class IceEnthalpySurface : public Diag<IceModel>
{
public:
  IceEnthalpySurface(const IceModel *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

//! \brief Computes enthalpy at the base of the ice.
class IceEnthalpyBasal : public Diag<IceModel>
{
public:
  IceEnthalpyBasal(const IceModel *m);
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
class IceHardness : public Diag<IceModel>
{
public:
  IceHardness(const IceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Effective viscosity of ice (3D). */
class IceViscosity : public Diag<IceModel>
{
public:
  IceViscosity(IceModel *m);
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
