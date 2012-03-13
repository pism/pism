// Copyright (C) 2012 David Maxwell
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

#ifndef _FEEVALUATOR_H_
#define _FEEVALUATOR_H_

#include "FETools.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"


class FEEvaluator {
public:
  FEEvaluator() : m_i(0), m_j(0),
  m_elementOwner(false), m_grid(NULL) { };
  virtual ~FEEvaluator() {};

  virtual void setElement( PetscInt i, PetscInt j);
  virtual void coordsAt( PetscInt q, PetscReal *x, PetscReal *y );
  virtual void testFunctionValues(PetscInt q, PetscInt k, PetscReal *v, PetscReal *dx, PetscReal *dy);
  virtual PetscReal JxW( PetscInt q );

protected:
  PetscInt  m_i, m_j;
  bool      m_elementOwner;
  IceGrid  *m_grid;
  
  FEQuadrature m_quadrature;
  FEDOFMap m_dofmap;
};

class FEEvaluator2S : public FEEvaluator {
public:
  FEEvaluator2S() : m_array(NULL), m_vec(NULL){ };
  virtual ~FEEvaluator2S();

  virtual PetscErrorCode init( IceModelVec2S &vec );
  virtual void release();
  
  virtual void setElement( PetscInt i, PetscInt j);
  virtual void valueAt( PetscInt q, PetscReal *v, PetscReal *dfdx, PetscReal *dfdy);
  virtual PetscReal nodeValue(PetscInt k);

protected:

  PetscReal m_values[FEQuadrature::Nq];
  PetscReal m_dfdx[FEQuadrature::Nq];
  PetscReal m_dfdy[FEQuadrature::Nq];
  
  PetscReal **m_array;
  IceModelVec2S *m_vec;
};  

class FEEvaluator2V : public FEEvaluator{
public:
  FEEvaluator2V() : m_array(NULL), m_vec(NULL) { };
  virtual ~FEEvaluator2V();

  virtual PetscErrorCode init( IceModelVec2V &vec );
  virtual void release();

  virtual void setElement( PetscInt i, PetscInt j);
  virtual void valueAt( PetscInt q, PetscReal *u, PetscReal *v, PetscReal *Du0, 
   PetscReal *Du1, PetscReal *Du2 );
  virtual PISMVector2 nodeValue(PetscInt k);


protected:

  PISMVector2  m_values[FEQuadrature::Nq];
  PetscReal m_Du[FEQuadrature::Nq][3];
  PISMVector2 **m_array;
  IceModelVec2V *m_vec;
};

#endif
