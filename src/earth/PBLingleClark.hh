/* Copyright (C) 2013, 2014 PISM Authors
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

#ifndef _PBLINGLECLARK_H_
#define _PBLINGLECLARK_H_

#include "PISMBedDef.hh"
#include <fftw3.h>
#include "deformation.hh"

namespace pism {

//! A wrapper class around BedDeformLC.
class PBLingleClark : public BedDef {
public:
  PBLingleClark(IceGrid &g, const Config &conf);
  virtual ~PBLingleClark();

  PetscErrorCode init(Vars &vars);
  PetscErrorCode update(double my_t, double my_dt);
protected:
  PetscErrorCode correct_topg();
  PetscErrorCode allocate();
  PetscErrorCode deallocate();
  PetscErrorCode transfer_to_proc0(IceModelVec2S *source, Vec result);
  PetscErrorCode transfer_from_proc0(Vec source, IceModelVec2S *result);
  Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
  VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
  // Vecs on processor 0:
  Vec Hp0,                      //!< ice thickness
    bedp0,                      //!< bed elevation
    Hstartp0,                   //!< initial (start-of-the-run) thickness
    bedstartp0,                 //!< initial bed elevation
    upliftp0;                   //!< bed uplift
  BedDeformLC bdLC;
};

} // end of namespace pism

#endif /* _PBLINGLECLARK_H_ */
