// Copyright (C) 2004--2015 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SIA_SLIDING_H_
#define _SIA_SLIDING_H_

#include "base/stressbalance/ShallowStressBalance.hh"

namespace pism {
namespace stressbalance {

/*!
 * This class implements an SIA sliding law.
 *
 * It is used by pismv test E \b only, hence the code duplication (the surface
 * gradient code is from SIAFD).
 */
class SIA_Sliding : public ShallowStressBalance {
public:
  SIA_Sliding(const IceGrid &g, const EnthalpyConverter &e);

  virtual ~SIA_Sliding();

  virtual void update(bool fast, const IceModelVec2S &melange_back_pressure);

protected:
  virtual void init_impl();

  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword,
                                       std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                     IO_Type nctype);

  virtual void compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual void surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual void surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual void surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual double basalVelocitySIA(double /*x*/, double /*y*/,
                                  double H, double T,
                                  double /*alpha*/, double mu,
                                  double min_T) const;
  IceModelVec2S m_work_2d;
  IceModelVec2Stag m_work_2d_stag[2]; // for the surface gradient
  double m_standard_gravity;

  bool m_verification_mode;
  std::string m_eisII_experiment;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SIA_SLIDING_H_ */
