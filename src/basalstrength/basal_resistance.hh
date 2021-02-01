// Copyright (C) 2004-2015, 2017, 2019 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#ifndef __basal_resistance_hh
#define __basal_resistance_hh

#include "pism/util/Units.hh"

namespace pism {

class Config;
class Logger;

//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
  This *pseudo* -plastic type can actually describe anything from linearly
  viscous till to purely plastic till.
*/
class IceBasalResistancePlasticLaw {
public:
  IceBasalResistancePlasticLaw(const Config &config);
  virtual ~IceBasalResistancePlasticLaw();
  virtual void print_info(const Logger &log, int threshold, units::System::Ptr system) const;
  virtual double drag(double tauc, double vx, double vy) const;
  virtual void drag_with_derivative(double tauc, double vx, double vy,
                                    double *drag, double *ddrag) const;
protected:
  double m_plastic_regularize;
};

class IceBasalResistancePseudoPlasticLaw : public IceBasalResistancePlasticLaw{
public:
  IceBasalResistancePseudoPlasticLaw(const Config &config);
  virtual ~IceBasalResistancePseudoPlasticLaw();
  virtual void print_info(const Logger &log, int threshold, units::System::Ptr system) const;
  virtual double drag(double tauc, double vx, double vy) const;
  virtual void drag_with_derivative(double tauc, double vx, double vy,
                                    double *drag, double *ddrag) const;
protected:
  double m_pseudo_q, m_pseudo_u_threshold, m_sliding_scale_factor_reduces_tauc;
};

class IceBasalResistanceRegularizedLaw : public IceBasalResistancePlasticLaw{
public:
  IceBasalResistanceRegularizedLaw(const Config &config);
  virtual ~IceBasalResistanceRegularizedLaw();
  virtual void print_info(const Logger &log, int threshold, units::System::Ptr system) const;
  virtual double drag(double tauc, double vx, double vy) const;
  virtual void drag_with_derivative(double tauc, double vx, double vy,
                                    double *drag, double *ddrag) const;
protected:
  double m_pseudo_q, m_pseudo_u_threshold, m_sliding_scale_factor_reduces_tauc;
};

} // end of namespace pism

#endif /* __basal_resistance_hh */

