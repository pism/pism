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


#ifndef SSATAUCTIKHONOVPROBLEM_HH_HB8UWICX
#define SSATAUCTIKHONOVPROBLEM_HH_HB8UWICX

#include <tr1/memory>

#include "TaoTikhonovProblem.hh"
#include "InvSSAForwardProblem.hh"


#include "TaoUtil.hh"
#include "functional/Functional.hh"


class SSATaucTikhonovProblem: public TaoTikhonovProblem<InvSSAForwardProblem> {
public:

  SSATaucTikhonovProblem( InvSSAForwardProblem &forward, 
                          SSATaucTikhonovProblem::DesignVec &d0, 
                          SSATaucTikhonovProblem::StateVec &u_obs, PetscReal eta, 
                          Functional<SSATaucTikhonovProblem::DesignVec>&designFunctional, 
                          Functional<SSATaucTikhonovProblem::StateVec>&stateFunctional) :
                TaoTikhonovProblem(forward,d0,u_obs,eta,designFunctional,stateFunctional) {};

  virtual ~SSATaucTikhonovProblem() {};

  virtual PetscErrorCode connect(TaoSolver tao);

  virtual PetscErrorCode getVariableBounds(TaoSolver tao, Vec lo, Vec hi); 

};

#endif /* end of include guard: SSATAUCTIKHONOVPROBLEM_HH_HB8UWICX */

