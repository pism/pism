// Copyright (C) 2012-2019, 2021, 2022 PISM Authors
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
  Distributed(std::shared_ptr<const IceGrid> g);
  virtual ~Distributed() = default;

  const array::Scalar& subglacial_water_pressure() const;

protected:
  virtual void restart_impl(const File &input_file, int record);

  virtual void bootstrap_impl(const File &input_file,
                              const array::Scalar &ice_thickness);

  virtual void init_impl(const array::Scalar &W_till,
                               const array::Scalar &W,
                               const array::Scalar &P);

  virtual double max_timestep_P_diff(double phi0, double dt_diff_w) const;

  void update_impl(double t, double dt, const Inputs& inputs);

  std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  void check_P_bounds(array::Scalar &P,
                      const array::Scalar &P_o,
                      bool enforce_upper);

  void P_from_W_steady(const array::Scalar &W,
                       const array::Scalar &P_overburden,
                       const array::Scalar &sliding_speed,
                       array::Scalar &result);

  void update_P(double dt,
                const array::CellType &cell_type,
                const array::Scalar &sliding_speed,
                const array::Scalar &surface_input_rate,
                const array::Scalar &basal_melt_rate,
                const array::Scalar &P_overburden,
                const array::Scalar &Wtill,
                const array::Scalar &Wtill_new,
                const array::Scalar &P,
                const array::Scalar1 &W,
                const array::Staggered1 &Ws,
                const array::Staggered1 &K,
                const array::Staggered1 &Q,
                array::Scalar &P_new) const;
protected:
  array::Scalar1 m_P;
  array::Scalar m_Pnew;
private:
  void initialization_message() const;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* _DISTRIBUTED_H_ */
