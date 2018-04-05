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

#ifndef _DISTRIBUTED_H_
#define _DISTRIBUTED_H_

#include "Routing.hh"

namespace pism {

namespace hydrology {

//! \brief The PISM subglacial hydrology model for a distributed linked-cavity system.
/*!
  This class implements the model documented in [\ref BuelervanPeltDRAFT].

  Unlike hydrology::Routing, the water pressure \f$P\f$ is a state variable, and there
  are modeled mechanisms for cavity geometry evolution, including creep closure
  and opening through sliding ("cavitation").  Because of cavitation, this model
  needs access to a StressBalance object.   Background references for this kind of
  model includes especially [\ref Kamb1987, \ref Schoofetal2012], but see also
  [\ref Hewitt2011, \ref Hewittetal2012, \ref Hewitt2013].

  In addition to the actions within the null strip taken by hydrology::Routing,
  this model also sets the staggered grid values of the gradient of the hydraulic
  potential to zero if either regular grid neighbor is in the null strip.
*/
class Distributed : public Routing {
public:
  Distributed(IceGrid::ConstPtr g);
  virtual ~Distributed();

  const IceModelVec2S& subglacial_water_pressure() const;

protected:
  virtual void restart_impl(const PIO &input_file, int record);

  virtual void bootstrap_impl(const PIO &input_file,
                              const IceModelVec2S &ice_thickness);

  virtual void initialize_impl(const IceModelVec2S &W_till,
                               const IceModelVec2S &W,
                               const IceModelVec2S &P);

  virtual double max_timestep_P_diff(double phi0, double dt_diff_w) const;

  void update_impl(double t, double dt, const Inputs& inputs);

  std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  void check_P_bounds(IceModelVec2S &P,
                      const IceModelVec2S &P_o,
                      bool enforce_upper);

  void P_from_W_steady(const IceModelVec2S &W,
                       const IceModelVec2S &P_overburden,
                       const IceModelVec2S &sliding_speed,
                       IceModelVec2S &result);

  void update_P(double dt,
                const IceModelVec2CellType &cell_type,
                const IceModelVec2S &sliding_speed,
                const IceModelVec2S &total_input,
                const IceModelVec2S &P_overburden,
                const IceModelVec2S &Wtill,
                const IceModelVec2S &Wtill_new,
                const IceModelVec2S &P,
                const IceModelVec2S &W,
                const IceModelVec2Stag &Ws,
                const IceModelVec2Stag &K,
                const IceModelVec2Stag &Q,
                IceModelVec2S &P_new) const;
protected:
  IceModelVec2S m_P;
  IceModelVec2S m_Pnew;
private:
  void initialization_message() const;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* _DISTRIBUTED_H_ */
