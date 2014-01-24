// Copyright (C) 2004-2014 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#include <petscsys.h>

#include "PISMUnits.hh"

class PISMConfig;

//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
This \e pseudo -plastic type can actually describe anything from linearly 
viscous till to purely plastic till.
 */
class IceBasalResistancePlasticLaw {
public:
  IceBasalResistancePlasticLaw(const PISMConfig &config);
  virtual ~IceBasalResistancePlasticLaw() {}
  virtual PetscErrorCode print_info(int verbthresh, MPI_Comm com);
  virtual double drag(double tauc, double vx, double vy);
  //! The derivative of drag with respect to \f$ alpha=\frac 1 2 |u|^2 \f$.
  virtual void drag_with_derivative(double tauc, double vx, double vy,
                                    double *drag, double *ddrag) const;
protected:
  double plastic_regularize;
  PISMUnitSystem m_unit_system;
};

class IceBasalResistancePseudoPlasticLaw : public IceBasalResistancePlasticLaw{
public:
  IceBasalResistancePseudoPlasticLaw(const PISMConfig &config);
  virtual ~IceBasalResistancePseudoPlasticLaw() {}
  virtual PetscErrorCode print_info(int verbthresh, MPI_Comm com);
  virtual double drag(double tauc, double vx, double vy);
  virtual void drag_with_derivative(double tauc, double vx, double vy,
                                    double *drag, double *ddrag) const;
protected:
  double pseudo_q, pseudo_u_threshold, sliding_scale;
};

#endif /* __basal_resistance_hh */

