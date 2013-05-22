// Copyright (C) 2012  David Maxwell
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

#ifndef IP_L2NORMFUNCTIONAL_HH_BSVF8BMQ
#define IP_L2NORMFUNCTIONAL_HH_BSVF8BMQ

#include "IPFunctional.hh"

class IP_L2NormFunctional2S : public IPInnerProductFunctional<IceModelVec2S> {
public:
  IP_L2NormFunctional2S(IceGrid &grid) : IPInnerProductFunctional<IceModelVec2S>(grid) {};
  virtual ~IP_L2NormFunctional2S() {};
  
  virtual PetscErrorCode valueAt(IceModelVec2S &x, PetscReal *OUTPUT);
  virtual PetscErrorCode dot(IceModelVec2S &a, IceModelVec2S &b, PetscReal *v);
  virtual PetscErrorCode gradientAt(IceModelVec2S &x, IceModelVec2S &gradient);

private:
  IP_L2NormFunctional2S(IP_L2NormFunctional2S const &);
  IP_L2NormFunctional2S & operator=(IP_L2NormFunctional2S const &);  
};

class IP_L2NormFunctional2V : public IPInnerProductFunctional<IceModelVec2V> {
public:
  IP_L2NormFunctional2V(IceGrid &grid) : IPInnerProductFunctional<IceModelVec2V>(grid) {};
  virtual ~IP_L2NormFunctional2V() {};
  
  virtual PetscErrorCode valueAt(IceModelVec2V &x, PetscReal *v);
  virtual PetscErrorCode dot(IceModelVec2V &a, IceModelVec2V &b, PetscReal *v);
  virtual PetscErrorCode gradientAt(IceModelVec2V &x, IceModelVec2V &gradient);

private:
  IP_L2NormFunctional2V(IP_L2NormFunctional2V const &);
  IP_L2NormFunctional2V & operator=(IP_L2NormFunctional2V const &);  
};

#endif /* end of include guard: IP_L2NORMFUNCTIONAL_HH_BSVF8BMQ */
