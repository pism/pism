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

#include "pism/coupler/FrontalMeltModel.hh"

namespace pism {
namespace frontalmelt {

/*!
 * Base class for frontal melt models with dedicated storage for output fields.
 *
 * All ocean models have storage for melange back pressure. (All but one set it to zero.)
 */
class CompleteFrontalMeltModel : public FrontalMeltModel {
public:
  // "modifier" constructor
  CompleteFrontalMeltModel(IceGrid::ConstPtr g, std::shared_ptr<FrontalMeltModel> input);
  // "model" constructor
  CompleteFrontalMeltModel(IceGrid::ConstPtr g);

  virtual ~CompleteFrontalMeltModel();
protected:
  virtual const IceModelVec2S& frontal_melt_rate_impl() const;

  IceModelVec2S::Ptr m_frontal_melt_rate;
};

} // end of namespace frontalmelt
} // end of namespace pism

#endif /* COMPLETEOCEANMODEL_H */
