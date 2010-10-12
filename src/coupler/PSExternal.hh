#include "PISMSurface.hh"

//! \brief A derived class created to couple PISM to an energy balance model (through
//! files).
/*!
  Uses an atmospheric lapse rate to correct temperatures at the ice surface.
 */ 
class PSExternal : public PISMSurfaceModel {
public:
  PSExternal(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  {
    gamma = 0;                  // essentially disables the lapse rate correction
    update_interval = 1;        // years
  }

  virtual ~PSExternal() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2S &result);
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
protected:
  PetscReal gamma, update_interval;
  IceModelVec2S acab, artm, artm_0, *usurf, *topg;
  string in_filename,           //!< the name of the file to read acab from
    out_filename;               //!< the name of the file to write

  virtual PetscErrorCode lock(string name, int op, int &fd);
  virtual PetscErrorCode unlock(int fd);
  virtual PetscErrorCode update_artm();
  virtual PetscErrorCode update_acab();
  virtual PetscErrorCode write_coupling_fields();
};

