/* Copyright (C) 2014 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef _PSVERIFICATION_H_
#define _PSVERIFICATION_H_

#include "PSFormulas.hh"

namespace pism {

class EnthalpyConverter;

//! Climate inputs for verification tests.
class PSVerification : public PSFormulas {
public:
  PSVerification(IceGrid &g, EnthalpyConverter *EC, int test);
  ~PSVerification();

  // the interface:
  void init();
  void update(PetscReal t, PetscReal dt);
private:
  int m_testname;
  EnthalpyConverter *m_EC;
  void update_ABCDEH(double t);
  void update_FG(double t);
  void update_KO();
  void update_L();
  void update_V();

  // FIXME: get rid of code duplication (these constants are defined
  // both here and in IceCompModel).
  static const double secpera, ablationRateOutside;

  static const PetscScalar ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static const PetscScalar Tmin;    // K;       minimum temperature (at center)
  static const PetscScalar LforFG;  // m;  exact radius of tests F&G ice sheet
  static const PetscScalar ApforG;  // m;  magnitude A_p of annular perturbation for test G;

};

} // end of namespace pism

#endif /* _PSVERIFICATION_H_ */
