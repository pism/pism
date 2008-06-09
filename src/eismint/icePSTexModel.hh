// Copyright (C) 2007-2008 Ed Bueler
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

#ifndef __icePSTexModel_hh
#define __icePSTexModel_hh

#include <petsc.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "iceEISModel.hh"

//! Derived class for a Plastic till ice Stream with Thermocoupling experiment.
/*!
This derived class supercedes an older class IceEISplModel.  The results from
that model were presented by Bueler at AGU 2007 and at NYU in Feb 2008.  The new 
model is more efficient in doing parameter studies, by having three streams going 
simultaneously.  On the other hand there are less command line options.  The
configuration is hard-wired into this derived class.
Also there is no "lake" or "fjord".
 */
class IcePSTexModel : public IceEISModel {

public:
  IcePSTexModel(IceGrid &g, IceType *i);
  virtual PetscErrorCode initFromOptions();
    
protected:
  char exper_chosen_name[10];
  
  PetscErrorCode setBedElev();
  PetscErrorCode setTillPhi();
  virtual PetscErrorCode summaryPrintLine(
              const PetscTruth printPrototype, const PetscTruth tempAndAge,
              const PetscScalar year, const PetscScalar dt, 
              const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
              const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);
private:
  int exper_chosen;
  int sectorNumberP2(const PetscScalar x, const PetscScalar y);
  bool inStream(const PetscScalar angle, const PetscScalar width,
               const PetscScalar x, const PetscScalar y,
               PetscScalar &x_loc, PetscScalar &y_loc);
  bool inStreamNbhd(bool strictly_in_stream,
               const PetscScalar angle, const PetscScalar width,
               const PetscScalar x, const PetscScalar y,
               PetscScalar &x_loc, PetscScalar &y_loc);
/*
  int inStreamP2(const PetscScalar width, 
                 const PetscScalar x, const PetscScalar y,
                 PetscScalar &x_loc, PetscScalar &y_loc);
  int inStreamNbhdP2(bool strictly_in_stream,
                 const PetscScalar width, 
                 const PetscScalar x, const PetscScalar y,
                 PetscScalar &x_loc, PetscScalar &y_loc);
*/
  PetscScalar phiLocal(const PetscScalar width, 
         const PetscScalar x, const PetscScalar y,
         const PetscScalar STRONG, const PetscScalar UP, const PetscScalar DOWN);
};

#endif /* __icePSTexModel_hh */

