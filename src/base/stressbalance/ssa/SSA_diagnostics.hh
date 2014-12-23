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

#ifndef _SSA_DIAGNOSTICS_H_
#define _SSA_DIAGNOSTICS_H_

#include "SSA.hh"

namespace pism {

//! \brief Computes the magnitude of the driving shear stress at the base of
//! ice (diagnostically).
class SSA_taud_mag : public Diag<SSA>
{
public:
  SSA_taud_mag(SSA *m);
  virtual void compute(IceModelVec* &result);
};

//! @brief Computes the driving shear stress at the base of ice
//! (diagnostically).
/*! This is *not* a duplicate of SSB_taud: SSA_taud::compute() uses
  SSA::compute_driving_stress(), which tries to be smarter near ice margins.
*/
class SSA_taud : public Diag<SSA>
{
public:
  SSA_taud(SSA *m);
  virtual void compute(IceModelVec* &result);
};


} // end of namespace pism

#endif /* _SSA_DIAGNOSTICS_H_ */
