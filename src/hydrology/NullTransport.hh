// Copyright (C) 2012-2018 PISM Authors
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

#ifndef _NULLTRANSPORT_H_
#define _NULLTRANSPORT_H_

#include "Hydrology.hh"

namespace pism {

namespace hydrology {

//! The PISM minimal model has till, but water that exceeds the capacity of the till is
//! not conserved. There is no model for lateral transport.
/*!
  This is the minimum functional derived class.  It updates till water thickness.
  It implements a version of the "undrained plastic bed" model of [\ref Tulaczyketal2000b],
  but with non-conserved drainage.

  This model can give no meaningful report on conservation errors.

  This talk illustrates a "till-can" metaphor applicable to this model:
  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf
*/
class NullTransport : public Hydrology {
public:
  NullTransport(IceGrid::ConstPtr g);
  virtual ~NullTransport();

protected:
  virtual void restart_impl(const PIO &input_file, int record);

  virtual void bootstrap_impl(const PIO &input_file,
                              const IceModelVec2S &ice_thickness);

  virtual void initialize_impl(const IceModelVec2S &W_till,
                               const IceModelVec2S &W,
                               const IceModelVec2S &P);

  virtual MaxTimestep max_timestep_impl(double t) const;

  //! Solves an implicit step of a highly-simplified ODE.
  virtual void update_impl(double t, double dt, const Inputs& inputs);

  void diffuse_till_water(double dt);

private:
  double m_diffuse_tillwat;
  double m_diffusion_time;
  double m_diffusion_distance;
  double m_tillwat_max;
  double m_tillwat_decay_rate;

  IceModelVec2S m_Wtill_old;

  void initialization_message() const;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* _NULLTRANSPORT_H_ */
