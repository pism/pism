#ifndef _IBSURFACEMODEL_HH_
#define _IBSURFACEMODEL_HH_

#include <coupler/PISMSurface.hh>
#include <base/util/iceModelVec.hh>
#include <coupler/PISMAtmosphere.hh>
#include <base/util/PISMVars.hh>
#include <coupler/surface/PSConstantPIK.hh>

namespace pism {
namespace icebin {

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-GLINT2, using a latitude
//! and surface elevation-dependent formula.

class IBSurfaceModel : public pism::surface::PIK {
public:

  /** @param conf Not Used (Looked up all the constructors, it just
      sets this->config, whic his not used
      @param g glint2::IceGrid*/
  IBSurfaceModel(IceGrid::ConstPtr grid);

  virtual ~IBSurfaceModel() {}

protected:
  void init_impl();
  void update_impl(double my_t, double my_dt);
};

}}		// namespace

#endif /* _IBSURFACEMODEL_HH_ */
