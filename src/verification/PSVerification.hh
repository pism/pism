/* Copyright (C) 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/coupler/surface/Formulas.hh"
#include "pism/util/EnthalpyConverter.hh"

namespace pism {
namespace surface {

//! Climate inputs for verification tests.
class Verification : public PSFormulas {
public:
  Verification(IceGrid::ConstPtr g, EnthalpyConverter::Ptr EC, int test);
  ~Verification();
private:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  MaxTimestep max_timestep_impl(double t) const;

  int m_testname;
  EnthalpyConverter::Ptr m_EC;
  void update_ABCDH(double t);
  void update_FG(double t);
  void update_KO();
  void update_L();
  void update_V();

  // FIXME: get rid of code duplication (these constants are defined
  // both here and in IceCompModel).
  static const double secpera, ablationRateOutside;

  static const double ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static const double Tmin;    // K;       minimum temperature (at center)
  static const double LforFG;  // m;  exact radius of tests F&G ice sheet
  static const double ApforG;  // m;  magnitude A_p of annular perturbation for test G;

};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSVERIFICATION_H_ */
