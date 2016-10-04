// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"

#include "base/util/IceGrid.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec.hh"
#include "base/age/AgeColumnSystem.hh"
#include "base/age/AgeModel.hh"

namespace pism {

//! Take a semi-implicit time-step for the age equation.
/*!
Let \f$\tau(t,x,y,z)\f$ be the age of the ice.  Denote the three-dimensional
velocity field within the ice fluid as \f$(u,v,w)\f$.  The age equation
is \f$d\tau/dt = 1\f$, that is, ice may move but it gets one year older in one
year.  Thus
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1 \f]
This equation is purely advective and hyperbolic.  The right-hand side is "1" as
long as age \f$\tau\f$ and time \f$t\f$ are measured in the same units.
Because the velocity field is incompressible, \f$\nabla \cdot (u,v,w) = 0\f$,
we can rewrite the equation as
    \f[ \frac{\partial \tau}{\partial t} + \nabla \left( (u,v,w) \tau \right) = 1 \f]
There is a conservative first-order numerical method; see AgeColumnSystem::solveThisColumn().

The boundary condition is that when the ice falls as snow it has age zero.
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas.  There is no
boundary condition elsewhere on the ice upper surface, as the characteristics
go outward in the ablation zone.  If the velocity in the bottom cell of ice
is upward (\f$w>0\f$) then we also apply a zero age boundary condition,
\f$\tau(t,x,y,0) = 0\f$.  This is the case where ice freezes on at the base,
either grounded basal ice freezing on stored water in till, or marine basal ice.
(Note that the water that is frozen-on as ice might be quite "old" in the sense
that its most recent time in the atmosphere was long ago; this comment is
relevant to any analysis which relates isotope ratios to modeled age.)

The numerical method is a conservative form of first-order upwinding, but the
vertical advection term is computed implicitly.  Thus there is no CFL-type
stability condition from the vertical velocity; CFL is only for the horizontal
velocity.  We use a finely-spaced, equally-spaced vertical grid in the
calculation.  Note that the columnSystemCtx methods coarse_to_fine() and
fine_to_coarse() interpolate back and forth between this fine grid and
the storage grid.  The storage grid may or may not be equally-spaced.  See
AgeColumnSystem::solveThisColumn() for the actual method.
 */
void IceModel::ageStep(const AgeModelInputs &inputs, double dt) {

  inputs.check();

  const IceModelVec2S &ice_thickness = *inputs.ice_thickness;

  const IceModelVec3
    &u3 = *inputs.u3,
    &v3 = *inputs.v3,
    &w3 = *inputs.w3;

  IceGrid::ConstPtr grid = m_ice_age.get_grid();

  AgeColumnSystem system(grid->z(), "age",
                         grid->dx(), grid->dy(), dt,
                         m_ice_age, u3, v3, w3); // linear system to solve in each column

  size_t Mz_fine = system.z().size();
  std::vector<double> x(Mz_fine);   // space for solution

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(m_ice_age);
  list.add(m_work3d);

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      system.init(i, j, ice_thickness(i, j));

      if (system.ks() == 0) {
        // if no ice, set the entire column to zero age
        m_work3d.set_column(i, j, 0.0);
      } else {
        // general case: solve advection PDE

        // solve the system for this column; call checks that params set
        system.solve(x);

        // put solution in IceModelVec3
        system.fine_to_coarse(x, i, j, m_work3d);

        // Ensure that the age of the ice is non-negative.
        //
        // FIXME: this is a kludge. We need to ensure that our numerical method has the maximum
        // principle instead. (We may still need this for correctness, though.)
        double *column = m_work3d.get_column(i, j);
        for (unsigned int k = 0; k < grid->Mz(); ++k) {
          if (column[k] < 0.0) {
            column[k] = 0.0;
          }
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_work3d.update_ghosts(m_ice_age);
}

} // end of namespace pism
