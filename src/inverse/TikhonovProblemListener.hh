// Copyright (C) 2012  David Maxwell
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

#ifndef TIKHONOVPROBLEMLISTENER_HH_WT5F3KD9
#define TIKHONOVPROBLEMLISTENER_HH_WT5F3KD9

#include <tr1/memory>

template<class InverseProblem> class TikhonovProblem;

template<class InverseProblem>
class TikhonovProblemListener {
public:
  typedef std::tr1::shared_ptr<TikhonovProblemListener<InverseProblem> > Ptr;
  typedef typename InverseProblem::DesignVec DesignVec;
  typedef typename InverseProblem::StateVec StateVec;
  
  TikhonovProblemListener() {}
  virtual ~TikhonovProblemListener() {}
  
  virtual PetscErrorCode 
  iteration( TikhonovProblem<InverseProblem> &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,  StateVec &diff_u,  DesignVec &grad_u,
             DesignVec &gradient) = 0;
};

template<class InverseProblem>
class SimpleTikhonovListener: public TikhonovProblemListener<InverseProblem> {
public:
  typedef typename InverseProblem::DesignVec DesignVec;
  typedef typename InverseProblem::StateVec StateVec;

  virtual ~SimpleTikhonovListener() {};
  
  virtual PetscErrorCode 
  iteration( TikhonovProblem<InverseProblem> &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,  StateVec &diff_u,  DesignVec &grad_u,
             DesignVec &gradient) {
    PetscErrorCode ierr;
    ierr = verbPrintf(1,PETSC_COMM_WORLD,"Tikhonov interation %d: objective %g penalty %g\n",iter,objectiveValue,designValue); CHKERRQ(ierr);
    return 0;
  }
  
};

#endif /* end of include guard: TIKHONOVPROBLEMLISTENER_HH_WT5F3KD9 */
