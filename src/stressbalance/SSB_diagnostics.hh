/* Copyright (C) 2014, 2015, 2016, 2017, 2023, 2026 PISM Authors
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

#ifndef PISM_SSB_DIAGNOSTICS_HH
#define PISM_SSB_DIAGNOSTICS_HH

#include "pism/util/Diagnostic.hh"

namespace pism {

class Vars;

namespace stressbalance {

class SSB_beta : public Diag<ShallowStressBalance>
{
public:
  SSB_beta(const ShallowStressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_SSB_DIAGNOSTICS_HH */
