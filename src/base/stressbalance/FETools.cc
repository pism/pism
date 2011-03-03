// Copyright (C) 2009--2011 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

// Utility functions used by the SSAFEM code.

#include "FETools.hh"

const FEShapeQ1::ShapeFunctionSpec FEShapeQ1::shapeFunction[FEShapeQ1::Nk] = 
{FEShapeQ1::shape0, FEShapeQ1::shape1, FEShapeQ1::shape2, FEShapeQ1::shape3};


FEElementMap::FEElementMap(const IceGrid &g)
{
  // Start by assuming ghost elements exist in all directions.
  // If there are 'mx' by 'my' points in the grid and no 
  // ghosts, there are (mx-1) by (my-1) elements.  In the presence
  // of ghosts there are then (mx+1) by (my+1) elements.
  xs= g.xs-1; xm = g.xm+1;
  ys= g.ys-1; ym = g.ym+1;

  // Now correct if needed. The only way there will not be ghosts is if the 
  // grid is not periodic and we are up against the grid boundary.
  
  if( !(g.periodicity & X_PERIODIC) )
  {
    // First element has x-index 0.
    if(xs < 0){
      xs = 0;
    }
    // Total number of elements is Mx-1, so xs+xm should be no larger.
    if(xs+xm > g.Mx-1) {
      xm = g.Mx-1-xs;
    }
  }

  if( !(g.periodicity & Y_PERIODIC) )
  {
    if(ys < 0){
      ys = 0;
    }
    if(ys+ym > g.My-1) {
      ym = g.My-1-ys;
    }
  }
}


/*! \brief Extract local degrees of freedom for element (\a i,\a j) from global vector \a xg to
  local vector \a x (scalar-valued DOF version). */
void FEDOFMap::extractLocalDOFs(PetscInt i,PetscInt j, PetscReal const*const*xg,PetscReal *x) const
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];  
}

/*! \brief Extract local degrees of freedom for element (\a i,\a j) from global vector \a xg to
local vector \a x (vector-valued DOF version).
*/
void FEDOFMap::extractLocalDOFs(PetscInt i,PetscInt j, PISMVector2 const*const*xg,PISMVector2 *x) const
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];    
}


//! Extract scalar degrees of freedom for the element specified previously with FEDOFMap::reset
void FEDOFMap::extractLocalDOFs(PetscReal const*const*xg,PetscReal *x) const
{
  extractLocalDOFs(m_i,m_j,xg,x);
}

//! Extract vector degrees of freedom for the element specified previously with FEDOFMap::reset
void FEDOFMap::extractLocalDOFs(PISMVector2 const*const*xg,PISMVector2 *x) const
{
  extractLocalDOFs(m_i,m_j,xg,x);
}

//! Convert a local degree of freedom index \a k to a global degree of freedom index (\a i,\a j).
void FEDOFMap::localToGlobal(PetscInt k, PetscInt *i, PetscInt *j)
{
  *i = m_i + kIOffset[k];
  *j = m_j + kJOffset[k];  
}

/*!\brief Initialize the FEDOFMap to element (\a i, \a j) for the purposes of inserting into
global residual and Jacobian arrays. */
void FEDOFMap::reset(PetscInt i, PetscInt j, const IceGrid &grid)
{
  m_i = i; m_j = j;
  // The meaning of i and j for a PISM IceGrid and for a Petsc DA are swapped (the so-called
  // fundamental transpose.  The interface between PISM and Petsc is the stencils, so all
  // interactions with the stencils involve a transpose.
  m_col[0].j = i;   m_col[0].i = j;
  m_col[1].j = i+1; m_col[1].i = j;
  m_col[2].j = i+1; m_col[2].i = j+1;
  m_col[3].j = i;   m_col[3].i = j+1;

  memcpy(m_row,m_col,Nk*sizeof(m_col[0]));

  // We do not ever sum into rows that are not owned by the local rank.
  for(PetscInt k=0; k<Nk; k++) {         
    PetscInt pism_i = m_row[k].j, pism_j = m_row[k].i;
    if (  pism_i < grid.xs || grid.xs+grid.xm-1 < pism_i || 
                pism_j < grid.ys || grid.ys+grid.ym-1 < pism_j ) {
      markRowInvalid(k);      
    }
  }
}

/*!\brief Mark that the row corresponding to local degree of freedom \a k should not be updated
when inserting into the global residual or Jacobian arrays. */
void FEDOFMap::markRowInvalid(PetscInt k)
{
  m_row[k].i=m_row[k].j = kDofInvalid;
}

/*!\brief Mark that the column corresponding to local degree of freedom \a k should not be updated
when inserting into the global Jacobian arrays. */
void FEDOFMap::markColInvalid(PetscInt k)
{
  m_col[k].i=m_col[k].j = kDofInvalid;
}

/*!\brief Add the values of element-local residual contributions \a y to the global residual
 vector \a yg. */
/*! The element-local residual should be an array of Nk values.*/
void FEDOFMap::addLocalResidualBlock(const PISMVector2 *y, PISMVector2 **yg)
{
  for (int k=0; k<Nk; k++) {
    if (m_row[k].i == kDofInvalid || m_row[k].j == kDofInvalid) continue;
    yg[m_row[k].j][m_row[k].i].u += y[k].u;
    yg[m_row[k].j][m_row[k].i].v += y[k].v;
  }
}

//! Add the contributions of an element-local Jacobian to the global Jacobian vector.
/*! The element-local Jacobian should be givnen as a row-major array of Nk*Nk values in the
scalar case or (2Nk)*(2Nk) values in the vector valued case. */
PetscErrorCode FEDOFMap::addLocalJacobianBlock(const PetscReal *K, Mat J)
{
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,Nk,m_row,Nk,m_col,K,ADD_VALUES);CHKERRQ(ierr);  
  return 0;
}

//! Set a diagonal entry for global degree of freedom (\a i ,\a j ) in a Jacobian matrix
/*! This is an unhappy hack for supporting Dirichlet constrained degrees of freedom.
In the scalar valued case, \a K should point to a single value, and in the vector case,
it should point to 4 (=2x2) values for the (2x2) block correspoinding to to u-u, u-v, v-u, and v-v
interactions at grid point (\a i, \aj).  Sheesh.*/
PetscErrorCode FEDOFMap::setJacobianDiag(PetscInt i, PetscInt j, const PetscReal*K, Mat J)
{
  MatStencil row;
  row.i=j; row.j=i;
  PetscErrorCode ierr = MatSetValuesBlockedStencil(J,1,&row,1,&row,K,INSERT_VALUES);CHKERRQ(ierr);  
  return 0;
}

const PetscInt FEDOFMap::kIOffset[4] = {0,1,1,0};
const PetscInt FEDOFMap::kJOffset[4] = {0,0,1,1};


FEQuadrature::FEQuadrature()
{  
}

//! Obtain the weights \f$w_q\f$ for quadrature.
void FEQuadrature::getWeightedJacobian(PetscReal *jxw)
{
  for(int q=0;q<Nq;q++)
  {
    jxw[q] = m_jacobianDet * quadWeights[q];
  } 
}

//! Obtain the weights \f$w_q\f$ for quadrature.
void FEQuadrature::init(const IceGrid &grid)
{
  // Since we use uniform cartesian coordinates, the Jacobian is constant and diagonal on every element.
  // Note that the reference element is \f$ [-1,1]^2 \f$ hence the extra factor of 1/2.
  PetscReal jacobian_x = 0.5*grid.dx;///ref.Length();
  PetscReal jacobian_y = 0.5*grid.dy;///ref.Length();
  m_jacobianDet = jacobian_x*jacobian_y;

  FEShapeQ1 shape;
  for(int q=0; q<Nq; q++){
    for(int k=0; k<Nk; k++){
      // PetscInt qk = q*Nk+k;
      shape.eval(k,quadPoints[q][0],quadPoints[q][1],&m_germs[q][k]);
      m_germs[q][k].dx /= jacobian_x;
      m_germs[q][k].dy /= jacobian_y;
    }
  }
}

//! Return the values at all quadrature points of all shape functions.
//* The return value is an Nq by Nk array of FEFunctionGerms. */
const FEFunctionGerm (*FEQuadrature::testFunctionValues())[FEQuadrature::Nq]
{
  return m_germs;
}

//! Return the values of all shape functions at quadrature point \a q
//* The return value is an array of Nk FEFunctionGerms. */
const FEFunctionGerm *FEQuadrature::testFunctionValues(PetscInt q)
{
  return m_germs[q];
}

//! Return the values at quadrature point \a q of shape function \a k.
const FEFunctionGerm *FEQuadrature::testFunctionValues(PetscInt q, PetscInt k)
{
  return m_germs[q] + k;
}


/*! \brief Compute the values at the quadrature ponits of a scalar-valued 
finite-element function with element-local  degrees of freedom \a x.*/
/*! There should be room for FEQuadrature::Nq values in the output vector \a vals. */
void FEQuadrature::computeTrialFunctionValues(const PetscReal *x, PetscReal *vals)
{
  for (int q=0; q<Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0;
    for (int k=0; k<Nk; k++) {
      vals[q] += test[k].val * x[k];
    }
  }
}

/*! \brief Compute the values and first derivatives at the quadrature 
points of a scalar-valued finite-element function with element-local  
degrees of freedom \a x.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors \a vals, \a dx, 
and \a dy. */
void FEQuadrature::computeTrialFunctionValues(const PetscReal *x, PetscReal *vals, PetscReal *dx, PetscReal *dy)
{
  for (int q=0; q<Nq; q++) {
    const FEFunctionGerm *test = m_germs[q];
    vals[q] = 0; dx[q] = 0; dy[q] = 0;
    for (int k=0; k<Nk; k++) {
      vals[q] += test[k].val * x[k];
      dx[q]   += test[k].dx * x[k];
      dy[q]   += test[k].dy * x[k];
    }
  }
}

/*! \brief Compute the values at the quadrature points on element (\a i,\a j) 
of a scalar-valued finite-element function with global degrees of freedom \a x.*/
/*! There should be room for FEQuadrature::Nq values in the output vector \a vals. */
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, 
                                             PetscReal const*const*xg, PetscReal *vals)
{
  dof.extractLocalDOFs(i,j,xg,m_tmpScalar);
  computeTrialFunctionValues(m_tmpScalar,vals);
}

/*! \brief Compute the values and first derivatives at the quadrature points 
on element (\a i,\a j)  of a scalar-valued finite-element function with global degrees of freedom \a x.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors \a vals, \a dx, and \a dy. */
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof, PetscReal const*const*xg, 
                                 PetscReal *vals, PetscReal *dx, PetscReal *dy)
{
  dof.extractLocalDOFs(i,j,xg,m_tmpScalar);
  computeTrialFunctionValues(m_tmpScalar,vals,dx,dy);  
}

/*! \brief Compute the values at the quadrature points of a vector-valued 
finite-element function with element-local degrees of freedom \a x.*/
/*! There should be room for FEQuadrature::Nq values in the output vector \a vals. */
void FEQuadrature::computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals)
{
  for (int q=0; q<Nq; q++) {
    vals[q].u = 0; vals[q].v = 0;
    const FEFunctionGerm *test = m_germs[q];
    for (int k=0; k<Nk; k++) {
      vals[q].u += test[k].val * x[k].u;
      vals[q].v += test[k].val * x[k].v;
    }
  }  

}

/*! \brief Compute the values and symmetric gradient at the quadrature points of a vector-valued 
finite-element function with element-local degrees of freedom \a x.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors \a vals and \a Dv. 
Each entry of \a Dv is an array of three numbers: 
\f[\left[\frac{du}{dx},\frac{dv}{dy},\frac{1}{2}\left(\frac{du}{dy}+\frac{dv}{dx}\right)\right]\f].
*/
void FEQuadrature::computeTrialFunctionValues( const PISMVector2 *x,  PISMVector2 *vals, PetscReal (*Dv)[3] )
{
  for (int q=0; q<Nq; q++) {
    vals[q].u = 0; vals[q].v = 0;
    PetscReal *Dvq = Dv[q];
    Dvq[0]=0; Dvq[1]=0; Dvq[2]=0;
    const FEFunctionGerm *test = m_germs[q];
    for(int k=0; k<Nk; k++) {
      vals[q].u += test[k].val * x[k].u;
      vals[q].v += test[k].val * x[k].v;
      Dvq[0] += test[k].dx * x[k].u;
      Dvq[1] += test[k].dy * x[k].v;
      Dvq[2] += 0.5*(test[k].dy*x[k].u + test[k].dx*x[k].v);
    }
  }  
}

/*! \brief Compute the values at the quadrature points of a vector-valued 
finite-element function on element (\a i,\a j) with global degrees of freedom \a xg.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors \a vals. */
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof,                                              
                                             PISMVector2 const*const*xg, PISMVector2 *vals )
{
  dof.extractLocalDOFs(i,j,xg,m_tmpVector);
  computeTrialFunctionValues(m_tmpVector,vals);
}

/*! \brief Compute the values and symmetric gradient at the quadrature points of a vector-valued 
finite-element function on element (\a i,\a j) with global degrees of freedom \a xg.*/
/*! There should be room for FEQuadrature::Nq values in the output vectors \a vals and \a Dv. 
Each entry of \a Dv is an array of three numbers: 
\f[\left[\frac{du}{dx},\frac{dv}{dy},\frac{1}{2}\left(\frac{du}{dy}+\frac{dv}{dx}\right)\right]\f].
*/
void FEQuadrature::computeTrialFunctionValues( PetscInt i, PetscInt j, const FEDOFMap &dof,                                              
                                             PISMVector2 const*const*xg, PISMVector2 *vals, PetscReal (*Dv)[3] )
{
  dof.extractLocalDOFs(i,j,xg,m_tmpVector);
  computeTrialFunctionValues(m_tmpVector,vals,Dv);
}

//! The quadrature points on the reference square \f$x,y=\pm 1/\sqrt{3}\f$.
const PetscReal FEQuadrature::quadPoints[FEQuadrature::Nq][2] = 
                                          {{ -0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573,  0.57735026918962573 },
                                           { -0.57735026918962573,  0.57735026918962573 }};

//! The weights w_i for gaussian quadrature on the reference element with these quadrature points
const PetscReal FEQuadrature::quadWeights[FEQuadrature::Nq]  = {1,1,1,1};



//! Legacy code that needs to vanish. \todo Make it go away.
PetscTruth Floating(const IceFlowLaw &ice, PetscScalar ocean_rho,
                           PetscReal H, PetscReal bed)
{
  return ice.rho*H + ocean_rho*bed < 0 ? PETSC_TRUE : PETSC_FALSE;
}


//! Legacy code that needs to vanish. \todo Make it go away.
int PismIntMask(PetscScalar maskvalue) {
  return static_cast<int>(floor(maskvalue + 0.5));
}

