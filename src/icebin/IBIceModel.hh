#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>

#include <base/util/IceGrid.hh>
#include <base/iceModel.hh>

#include <base/util/pism_options.hh>
#include <coupler/atmosphere/PAFactory.hh>
#include <coupler/ocean/POFactory.hh>
#include <coupler/surface/PSFactory.hh>

#include <base/util/PISMTime.hh>
// --------------------------------
#include <icebin/IBSurfaceModel.hh>

namespace pism {
namespace icebin {


/** This is the ICEBIN customized version of PISM's ::IceModel class.

See https://github.com/pism/pism/issues/219

Here's a short term solution, though: create a new class
IcebinIceModel derived from IceModel and re-implement
IceModel::allocate_couplers(). In it, set
IceModel::external_surface_model and IceModel::external_ocean_model as
you see fit (IceModel will not de-allocate a surface (ocean) model if
it is set to true) and allocate PSConstantICEBIN. You might also want
to add IBIceModel::get_surface_model() which returns
IceModel::surface to get access to PSConstantICEBIN from outside of
IBIceModel.
*/

class IBIceModel : public pism::IceModel {
public:

  // see iceModel.cc for implementation of constructor and destructor:
  IBIceModel(IceGrid::Ptr grid, Context::Ptr context);
  virtual ~IBIceModel(); // must be virtual merely because some members are virtual

  virtual void allocate_couplers();

  /** @return Our instance of PSConstantICEBIN */
  IBSurfaceModel *ib_surface_model() {
    return dynamic_cast<IBSurfaceModel *>(surface);
  }
};

}}
