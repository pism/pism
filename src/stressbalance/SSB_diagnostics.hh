/* Copyright (C) 2014, 2015, 2016, 2017, 2023 PISM Authors
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

#ifndef _SSB_DIAGNOSTICS_H_
#define _SSB_DIAGNOSTICS_H_

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

//! \brief Computes the gravitational driving stress (diagnostically).
class SSB_taud : public Diag<ShallowStressBalance>
{
public:
  SSB_taud(const ShallowStressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

//! \brief Computes the magnitude of the gravitational driving stress
//! (diagnostically).
class SSB_taud_mag : public Diag<ShallowStressBalance>
{
public:
  SSB_taud_mag(const ShallowStressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

//! @brief Computes the basal shear stress @f$ \tau_b @f$.
class SSB_taub : public Diag<ShallowStressBalance>
{
public:
  SSB_taub(const ShallowStressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

//! \brief Computes the magnitude of the basal shear stress
//! (diagnostically).
class SSB_taub_mag : public Diag<ShallowStressBalance>
{
public:
  SSB_taub_mag(const ShallowStressBalance *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl() const;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSB_DIAGNOSTICS_H_ */
