// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2023, 2026 Constantine Khroulev
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

#ifndef PISM_STRESSBALANCE_DIAGNOSTICS_HH
#define PISM_STRESSBALANCE_DIAGNOSTICS_HH

#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/Diagnostic.hh"

namespace pism {
namespace stressbalance {

//! \brief Computes basal frictional heating.
class PSB_bfrict : public Diag<StressBalance>
{
public:
  PSB_bfrict(const StressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

//! \brief Reports the volumetric strain heating (3D).
class PSB_strainheat : public Diag<StressBalance>
{
public:
  PSB_strainheat(const StressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_STRESSBALANCE_DIAGNOSTICS_HH */
