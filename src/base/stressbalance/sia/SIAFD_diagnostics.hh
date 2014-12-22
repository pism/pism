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

#ifndef _SIAFD_DIAGNOSTICS_H_
#define _SIAFD_DIAGNOSTICS_H_

#include "SIAFD.hh"

namespace pism {

//! \brief Computes the multiplier \f$\theta\f$ in Schoof's (2003) theory of the
//! effect of bed roughness on the diffusivity of the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFD_schoofs_theta : public Diag<SIAFD>
{
public:
  SIAFD_schoofs_theta(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

//! \brief Computes the smoothed bed elevation from Schoof's (2003) theory of the
//! effect of bed roughness on the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFD_topgsmooth : public Diag<SIAFD>
{
public:
  SIAFD_topgsmooth(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

//! \brief Computes the thickness relative to the smoothed bed elevation in
//! Schoof's (2003) theory of the effect of bed roughness on the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFD_thksmooth : public Diag<SIAFD>
{
public:
  SIAFD_thksmooth(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

//! \brief Compute diffusivity of the SIA flow.
class SIAFD_diffusivity : public Diag<SIAFD>
{
public:
  SIAFD_diffusivity(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

//! \brief Compute diffusivity of the SIA flow (on the staggered grid).
class SIAFD_diffusivity_staggered : public Diag<SIAFD>
{
public:
  SIAFD_diffusivity_staggered(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

//! \brief Reports the x-component of the ice surface gradient on the staggered
//! grid as computed by SIAFD.
class SIAFD_h_x : public Diag<SIAFD>
{
public:
  SIAFD_h_x(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

//! \brief Reports the y-component of the ice surface gradient on the staggered
//! grid as computed by SIAFD.
class SIAFD_h_y : public Diag<SIAFD>
{
public:
  SIAFD_h_y(SIAFD *m, IceGrid &g);
  virtual void compute(IceModelVec* &result);
};

} // end of namespace pism

#endif /* _SIAFD_DIAGNOSTICS_H_ */
