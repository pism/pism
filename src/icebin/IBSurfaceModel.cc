#include <icebin/IBSurfaceModel.hh>
#include <base/util/io/PIO.hh>
#include <base/util/PISMVars.hh>
#include <base/util/IceGrid.hh>

namespace pism {
namespace icebin {

///// Constant-in-time surface model for accumulation,
///// ice surface temperature parameterized as in PISM-GLINT2 dependent on latitude and surface elevation


IBSurfaceModel::IBSurfaceModel(IceGrid::ConstPtr grid)
  : pism::surface::PIK(grid)
{}

void IBSurfaceModel::init_impl() {
  // No atmosphere to initialize (aka pism::SurfaceModel)
  // No files to read (aka PSConstantPIK)
}


void IBSurfaceModel::update_impl(double my_t, double my_dt)
{
  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  const IceModelVec2S
    &usurf = *m_grid->variables().get_2d_scalar("surface_altitude"),
    &lat   = *m_grid->variables().get_2d_scalar("latitude");

#if 0
  // This code from PIK needs to go
  IceModelVec::AccessList list;
  list.add(m_ice_surface_temp);
  list.add(usurf);
  list.add(lat);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    m_ice_surface_temp(i,j) = 273.15 + 30 - 0.0075 * usurf(i,j) - 0.68775 * lat(i,j)*(-1.0);
  }
#endif

  printf("IBSurfaceModel::update(%f) dumping variables\n", my_t);
  m_ice_surface_temp.dump("ice_surface_temp.nc");
  m_climatic_mass_balance.dump("climatic_mass_balance.nc");
  printf("IBSurfaceModel::update(%f) done dumping variables\n", my_t);

}

} // namespace icebin
} // namespace pism
