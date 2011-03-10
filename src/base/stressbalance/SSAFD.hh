// Copyright (C) 2004--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SSAFD_H_
#define _SSAFD_H_

#include "SSA.hh"
#include <petscksp.h>

//! PISM's SSA solver: the finite difference implementation
class SSAFD : public SSA
{
public:
  SSAFD(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
        const NCConfigVariable &c) :
    SSA(g,b,i,e,c)
  {
    allocate_fd();  // can't be done by allocate() since constructor is not virtual
  }

  virtual ~SSAFD()
  {
    deallocate_fd();
  }

  virtual PetscErrorCode init(PISMVars &vars);

protected:
  virtual PetscErrorCode allocate_fd();

  virtual PetscErrorCode deallocate_fd();

  virtual PetscErrorCode solve();

  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result);

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               PetscReal epsilon);

  virtual PetscErrorCode compute_nuH_norm(PetscReal &norm,
                                          PetscReal &norm_change);

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs);

  virtual PetscErrorCode writeSSAsystemMatlab();

  virtual PetscErrorCode update_nuH_viewers();

  // objects used internally
  IceModelVec2Stag hardness, nuH, nuH_old;
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSARHS;

  bool view_nuh;
  PetscViewer nuh_viewer;
  PetscInt nuh_viewer_size;
};

//! \brief A derived class implementing PIK's calving frong boundary conditions.
class SSAFD_PIK : public SSAFD
{
public:
  SSAFD_PIK(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
             const NCConfigVariable &c)
    : SSAFD(g, b, i, e, c) {}
  virtual ~SSAFD_PIK() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);
  virtual PetscErrorCode assemble_rhs(Vec rhs);
};


//! Constructs a new SSAFD
SSA * SSAFDFactory(IceGrid &, IceBasalResistancePlasticLaw &, 
                  IceFlowLaw &, EnthalpyConverter &, const NCConfigVariable &);

// //! Constructs a new SSAFD_PIK
// SSA * SSAFD_PIKFactory(IceGrid &, IceBasalResistancePlasticLaw &, 
//                        IceFlowLaw &, EnthalpyConverter &, const NCConfigVariable &);

#endif /* _SSAFD_H_ */

