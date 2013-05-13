/* Copyright (C) 2013 PISM Authors
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

#ifndef _PISMICEBERGREMOVER_H_
#define _PISMICEBERGREMOVER_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"

/*! \brief PISM iceberg remover */
/*!
 * Identifies and removes free-floating icebergs, which cause
 * well-posedness problems for stress solvers.
 *
 * Icebergs are, in this context, floating regions that are _not_
 * attached, through a chain of positive thickness ice-filled cells,
 * to at least one grounded cell.
 *
 * They cause the SSA operator to have a nontrivial null space.
 *
 * They are observed to cause unrealistically large velocities that
 * may affect ice velocities elsewhere.
 *
 * This class uses a serial connected component labeling algorithm to
 * remove "icebergs".
 */
class PISMIcebergRemover : public PISMComponent
{
public:
  PISMIcebergRemover(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMIcebergRemover();

  virtual PetscErrorCode init(PISMVars &vars);
  PetscErrorCode update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness);
protected:
  PetscErrorCode allocate();
  PetscErrorCode deallocate();

  PetscErrorCode transfer_to_proc0();
  PetscErrorCode transfer_from_proc0();

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO& nc);
  
  Vec m_g2, m_g2natural;  //!< global Vecs used to transfer data to/from processor 0.
  VecScatter m_scatter; //!< VecScatter used to transfer data to/from processor 0.
  Vec m_mask_p0;
  DM m_da2;
};

#endif /* _PISMICEBERGREMOVER_H_ */
