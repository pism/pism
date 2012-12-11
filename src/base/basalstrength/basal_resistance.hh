// Copyright (C) 2004-2012 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#ifndef __basal_resistance_hh
#define __basal_resistance_hh

#include <petscsys.h>
#include "NCVariable.hh"

//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
This \e pseudo -plastic type can actually describe anything from linearly 
viscous till to purely plastic till.
 */
class IceBasalResistancePlasticLaw {
public:
  IceBasalResistancePlasticLaw(const NCConfigVariable &config);
  virtual ~IceBasalResistancePlasticLaw() {}
  virtual PetscErrorCode printInfo(int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
  // Also get the derivative of drag with respect to \f$ alpha=\frac 1 2 \abs{u}^2 \f$.
  virtual void dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
                                  PetscScalar *drag, PetscScalar *ddrag) const;
protected:
  double plastic_regularize;
};

class IceBasalResistancePseudoPlasticLaw : public IceBasalResistancePlasticLaw{
public:
  IceBasalResistancePseudoPlasticLaw(const NCConfigVariable &config);
  virtual ~IceBasalResistancePseudoPlasticLaw() {}
  virtual PetscErrorCode printInfo(int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
  // Also get the derivative of drag with respect to \f$ alpha=\frac 1 2 \abs{u}^2 \f$.
  virtual void dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy,
                                  PetscScalar *drag, PetscScalar *ddrag) const;
protected:
  PetscReal pseudo_q, pseudo_u_threshold, sliding_scale;
};

#endif /* __basal_resistance_hh */

