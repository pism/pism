/* Copyright (C) 2021, 2022, 2023 PISM Authors
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
#ifndef ICEBERGREMOVERFEM_H
#define ICEBERGREMOVERFEM_H

#include "pism/frontretreat/util/IcebergRemover.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace calving {

/*!
 * Iceberg remover used with FEM-based stress balance solvers.
 *
 * In the FEM context we need to identify patches of elements that are not connected to
 * grounded ice. (Two elements are considered "connected" if they share a boundary.)
 */
class IcebergRemoverFEM : public IcebergRemover {
public:
  IcebergRemoverFEM(std::shared_ptr<const Grid> g);
private:
  void update_impl(const array::Scalar &bc_mask,
                   array::CellType1 &cell_type,
                   array::Scalar &ice_thickness);
  array::Scalar m_mask;
};

} // end of namespace calving
} // end of namespace pism

#endif /* ICEBERGREMOVERFEM_H */
