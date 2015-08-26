// Copyright (C) 2004-2011, 2013 -- 2015 Florian Ziemen, Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petscsys.h>

#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/pism_options.hh"
#include "columnSystem.hh"
#include "iceModel.hh"


namespace pism {

//! Tridiagonal linear system for vertical column of tracer (pure advection) problem.
class tracerSystemCtx : public columnSystemCtx {
public:
  tracerSystemCtx(const std::vector<double>& storage_grid,
               const std::string &my_prefix,
               double dx, double dy, double dt,
               const IceModelVec3 &tracer,
               const IceModelVec3 &u3,
               const IceModelVec3 &v3,
               const IceModelVec3 &w3);

  void initThisColumn(int i, int j, double thickness);

  void solveThisColumn(std::vector<double> &x, const double upper_bound, const double lower_bound);
protected:
  const IceModelVec3 &m_tracer3;
  double m_nu;
  std::vector<double> m_A, m_A_n, m_A_e, m_A_s, m_A_w;
};


tracerSystemCtx::tracerSystemCtx(const std::vector<double>& storage_grid,
                           const std::string &my_prefix,
                           double dx, double dy, double dt,
                           const IceModelVec3 &tracer,
                           const IceModelVec3 &u3,
                           const IceModelVec3 &v3,
                           const IceModelVec3 &w3)
  : columnSystemCtx(storage_grid, my_prefix, dx, dy, dt, u3, v3, w3),
    m_tracer3(tracer) {

  size_t Mz = m_z.size();
  m_A.resize(Mz);
  m_A_n.resize(Mz);
  m_A_e.resize(Mz);
  m_A_s.resize(Mz);
  m_A_w.resize(Mz);

  m_nu = m_dt / m_dz; // derived constant
}

void tracerSystemCtx::initThisColumn(int i, int j, double thickness) {
  init_column(i, j, thickness);
  if (m_ks == 0) {
    return;
  }

  coarse_to_fine(m_u3, i, j, &m_u[0]);
  coarse_to_fine(m_v3, i, j, &m_v[0]);
  coarse_to_fine(m_w3, i, j, &m_w[0]);

  coarse_to_fine(m_tracer3, m_i, m_j,   &m_A[0]);
  coarse_to_fine(m_tracer3, m_i, m_j+1, &m_A_n[0]);
  coarse_to_fine(m_tracer3, m_i+1, m_j, &m_A_e[0]);
  coarse_to_fine(m_tracer3, m_i, m_j-1, &m_A_s[0]);
  coarse_to_fine(m_tracer3, m_i-1, m_j, &m_A_w[0]);
}

//! Conservative first-order upwind scheme with implicit in the vertical: one column solve. -- see age...
/*!
Surface and basal boundary conditions are provided. No source term.
Most of the documentatino can be found in ageSystemCtx::solveThisColum
 */
	void tracerSystemCtx::solveThisColumn(std::vector<double> &x, const double upper_bound, const double lower_bound ) {

  TridiagonalSystem &S = *m_solver;
  // set up system: 0 <= k < m_ks
  for (unsigned int k = 0; k < m_ks; k++) {
    // do lowest-order upwinding, explicitly for horizontal
    S.RHS(k) =  (m_u[k] < 0 ?
                 m_u[k] * (m_A_e[k] -  m_A[k]) / m_dx :
                 m_u[k] * (m_A[k]  - m_A_w[k]) / m_dx);
    S.RHS(k) += (m_v[k] < 0 ?
                 m_v[k] * (m_A_n[k] -  m_A[k]) / m_dy :
                 m_v[k] * (m_A[k]  - m_A_s[k]) / m_dy);
    // note it is the age eqn: dage/dt = 1.0 and we have moved the hor.
    //   advection terms over to right:
    S.RHS(k) = m_A[k] - m_dt * (S.RHS(k));

    // do lowest-order upwinding, *implicitly* for vertical
    double AA = m_nu * m_w[k];
    if (k > 0) {
      if (AA >= 0) { // upward velocity
        S.L(k) = - AA;
        S.D(k) = 1.0 + AA;
        S.U(k) = 0.0;
      } else { // downward velocity; note  -AA >= 0
        S.L(k) = 0.0;
        S.D(k) = 1.0 - AA;
        S.U(k) = + AA;
      }
    } else { // k == 0 case
      // note L[0] is not used
      if (AA > 0) { // if strictly upward velocity apply boundary condition:
                    // age = 0 because ice is being added to base
        S.D(0) = 1.0;
        S.U(0) = 0.0;
        S.RHS(0) = lower_bound;
      } else { // downward velocity; note  -AA >= 0
        S.D(0) = 1.0 - AA;
        S.U(0) = + AA;
        // keep rhs[0] as is
      }
    }
  }  // done "set up system: 0 <= k < m_ks"

  // surface b.c. at m_ks
  if (m_ks > 0) {
    S.L(m_ks) = 0;
    S.D(m_ks) = 1.0;   // ignore U[m_ks]
    S.RHS(m_ks) = upper_bound;  // age zero at surface
  }

  // solve it
  try {
    S.solve(m_ks + 1, x);
  }
  catch (RuntimeError &e) {
    e.add_context("solving the tri-diagonal system (ageSystemCtx) at (%d, %d)\n"
                  "saving system to m-file... ", m_i, m_j);
    reportColumnZeroPivotErrorMFile(m_ks + 1);
    throw;
  }
  // x[k] contains age for k=0,...,ks, but set age of ice above (and
  // at) surface to zero years
  for (unsigned int k = m_ks + 1; k < x.size(); k++) {
    x[k] = upper_bound; // SET TO UPPER BOUND?
  }
}


//! Take a semi-implicit time-step for the tracer equation.
/*!
See ageStep()

 */
void IceModel::tracerStep() {
  PetscErrorCode  ierr;

  bool viewOneColumn = options::Bool("-view_sys",
                                     "save column system information to file");

  const IceModelVec3
    &u3 = stress_balance->velocity_u(),
    &v3 = stress_balance->velocity_v(),
    &w3 = stress_balance->velocity_w();

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(ice_surface_elevation);
  list.add(tracer_x); // HAS TO BE DONE FOR EVERY TRACER!
  list.add(tracer_y); // HAS TO BE DONE FOR EVERY TRACER!
  list.add(tracer_z); // HAS TO BE DONE FOR EVERY TRACER!
  list.add(tracer_t); // HAS TO BE DONE FOR EVERY TRACER!
  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(vWork3d);

  IceModelVec3* tracers[4];
	tracers[0]=&tracer_x;
	tracers[1]=&tracer_y;
	tracers[2]=&tracer_z;
	tracers[3]=&tracer_t;

	for (int t = 0; t<4; t++){
		{
		tracerSystemCtx system(m_grid->z(), "tracer",
													 m_grid->dx(), m_grid->dy(), dt_TempAge,
													 *tracers[t], u3, v3, w3); // linear system to solve in each column

		size_t Mz_fine = system.z().size();
		std::vector<double> x(Mz_fine);   // space for solution

		ParallelSection loop(m_grid->com);
		try {
			for (Points p(*m_grid); p; p.next()) {
				const int i = p.i(), j = p.j();

				system.initThisColumn(i, j, ice_thickness(i, j));

				if (system.ks() == 0) {
					// NO ICE -- WHAT DO YOU WANT TO DO? CURRENTLY SETTING EVERYTHING TO -9e9
					switch (t) {
					case 0:
						vWork3d.set_column(i, j, m_grid->x(i));
						break ;
					case 1:
						vWork3d.set_column(i, j, m_grid->y(j));
						break ;
					case 2:
						vWork3d.set_column(i, j, ice_surface_elevation(i,j));
						break ;
					case 3:
						vWork3d.set_column(i, j, m_time->current());
					}
				} else {
					// general case: solve advection PDE

					// solve the system for this column; call checks that params set
					switch (t) {
					case 0:
						system.solveThisColumn(x, m_grid->x(i), m_grid->x(i));
						break ;
					case 1:
						system.solveThisColumn(x, m_grid->y(j), m_grid->y(j));
						break ;
					case 2:
						system.solveThisColumn(x, ice_surface_elevation(i,j), ice_surface_elevation(i,j)-ice_thickness(i,j));
						break ;
					case 3:
						system.solveThisColumn(x, m_time->current(), m_time->current());
					}

					if (viewOneColumn && (i == id && j == jd)) {
						ierr = PetscPrintf(PETSC_COMM_SELF,
															 "\n"
															 "in tracerStep(): saving tracerSystemCtx at (i,j)=(%d,%d) to m-file... \n",
															 i, j);
						PISM_CHK(ierr, "PetscPrintf");

						system.viewColumnInfoMFile(x);
					}

					// put solution in IceModelVec3
					system.fine_to_coarse(x, i, j, vWork3d);

					// Ensure that the age of the ice is non-negative.
					//
					// FIXME: this is a kludge. We need to ensure that our numerical method has the maximum
					// principle instead. (We may still need this for correctness, though.)
					// double *column = vWork3d.get_column(i, j);
					// for (unsigned int k = 0; k < m_grid->Mz(); ++k) {
					// 	if (column[k] < 0.0) {
					// 		column[k] = 0.0;
					// 	}
					// }
				}
			}
		} catch (...) {
			loop.failed();
		}
		loop.check();

		vWork3d.update_ghosts(*tracers[t]); // spread the news
	}
}

}
} // end of namespace pism
