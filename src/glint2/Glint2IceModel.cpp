#include "Glint2IceModel.hpp"

namespace glint2 {
namespace pism {

Glint2IceModel::Glint2IceModel(IceGrid &g, PISMConfig &config, PISMConfig &overrides) :
	::IceModel(g, config, overrides)
{}

Glint2IceModel::~Glint2IceModel() {} // must be virtual merely because some members are virtual




PetscErrorCode Glint2IceModel::allocate_couplers()
{
	PetscErrorCode ierr;
	// Initialize boundary models:
	PAFactory pa(grid, config);
	PSFactory ps(grid, config);
	POFactory po(grid, config);
	PISMAtmosphereModel *atmosphere;

	ierr = PetscOptionsBegin(grid.com, "", "Options choosing PISM boundary models", "");
        PISM_PETSC_CHK(ierr, "PetscOptionsBegin");

#if 1
	// GLINT2-modified version
	if (surface == NULL) {
		surface = new PSConstantGLINT2(grid, config);
		external_surface_model = false;

		pa.create(atmosphere);
		surface->attach_atmosphere_model(atmosphere);
	}
#else
	// Original Version
	if (surface == NULL) {
		ps.create(surface);
		external_surface_model = false;

		pa.create(atmosphere);
		surface->attach_atmosphere_model(atmosphere);
	}
#endif

	if (ocean == NULL) {
		po.create(ocean);
		external_ocean_model = false;
	}
	ierr = PetscOptionsEnd();
        PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

	return 0;
}
}}
