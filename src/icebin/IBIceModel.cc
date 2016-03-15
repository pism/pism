#include <icebin/IBIceModel.hh>

namespace pism {
namespace icebin {

IBIceModel::IBIceModel(IceGrid::Ptr grid, Context::Ptr context) :
    IceModel(grid, context) {}

IBIceModel::~IBIceModel() {} // must be virtual merely because some members are virtual

void IBIceModel::allocate_couplers()
{
  // Initialize boundary models:
  atmosphere::Factory pa(m_grid);
  surface::Factory ps(m_grid);
  ocean::Factory po(m_grid);
  atmosphere::AtmosphereModel *atmosphere;

  // http://stackoverflow.com/questions/3825668/c-c-checking-for-null-pointer
  if (!surface) {
    m_log->message(2,
             "# Allocating a surface process model or coupler...\n");

    // surface = ps.create();
    surface = new IBSurfaceModel(grid());
    external_surface_model = false;

    atmosphere = pa.create();
    surface->attach_atmosphere_model(atmosphere);
  }

  if (!ocean) {
    m_log->message(2,
             "# Allocating an ocean model or coupler...\n");

    ocean = po.create();
    external_ocean_model = false;
  }
}

}}
