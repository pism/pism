// Copyright (C) 2013  David Maxwell
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

#ifndef IPGROUNDEDICEH1NORMFUNCTIONAL_HH_Q4IZKJOR
#define IPGROUNDEDICEH1NORMFUNCTIONAL_HH_Q4IZKJOR


#include "IPFunctional.hh"
#include "Mask.hh"

class IPGroundedIceH1NormFunctional2S : public IPInnerProductFunctional<IceModelVec2S> {
public:
  IPGroundedIceH1NormFunctional2S(IceGrid &grid, PetscReal cL2, 
      PetscReal cH1, IceModelVec2Int &ice_mask, IceModelVec2Int *dirichletLocations=NULL) :
      IPInnerProductFunctional<IceModelVec2S>(grid),
      m_cL2(cL2), m_cH1(cH1), m_dirichletIndices(dirichletLocations),  m_ice_mask(ice_mask) {};
  virtual ~IPGroundedIceH1NormFunctional2S() {};
  
  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

protected:

  PetscReal m_cL2, m_cH1;
  IceModelVec2Int *m_dirichletIndices;
  IceModelVec2Int &m_ice_mask;

private:
  IPGroundedIceH1NormFunctional2S(IPGroundedIceH1NormFunctional2S const &);
  IPGroundedIceH1NormFunctional2S & operator=(IPGroundedIceH1NormFunctional2S const &);  
};

#endif /* end of include guard: IPGROUNDEDICEH1NORMFUNCTIONAL_HH_Q4IZKJOR */
