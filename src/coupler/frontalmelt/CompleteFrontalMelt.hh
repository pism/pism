/* Copyright (C) 2018 PISM Authors
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

#ifndef COMPLETEFRONTALMELTMODEL_H
#define COMPLETEFRONTALMELTMODEL_H

#include "pism/coupler/FrontalMelt.hh"

namespace pism {
namespace frontalmelt {

/*!
 * Base class for frontal melt models with dedicated storage for output fields.
 *
 */
class CompleteFrontalMelt : public FrontalMelt {
public:
  // "modifier" constructor
  CompleteFrontalMelt(IceGrid::ConstPtr g, std::shared_ptr<FrontalMelt> input);
  // "model" constructor
  CompleteFrontalMelt(IceGrid::ConstPtr g);

  virtual ~CompleteFrontalMelt();
protected:
  virtual const IceModelVec2S& frontal_melt_rate_impl() const;

  IceModelVec2S::Ptr m_frontal_melt_rate;
};

} // end of namespace frontalmelt
} // end of namespace pism

#endif /* COMPLETEOCEANMODEL_H */
