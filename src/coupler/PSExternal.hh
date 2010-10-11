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

  virtual PetscErrorCode init(PISMVars &vars); // done
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }             // done

  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2S &result); // done
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2S &result); // done
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);           // needs work
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years); // done
protected:
  PetscReal gamma, update_interval;
  IceModelVec2S acab, artm, usurf_0, *usurf, *topg;
  string in_filename,           //!< the name of the file to read acab from
    out_filename;               //!< the name of the file to write

  PetscErrorCode lock(string name, int op, int &fd);
  PetscErrorCode unlock(int fd);
};

