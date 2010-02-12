#ifndef __PA_EISMINT_Greenland
#define __PA_EISMINT_Greenland

#include "../coupler/PISMAtmosphere.hh"

class PA_EISMINT_Greenland : public PAFausto {
public:
  PA_EISMINT_Greenland(IceGrid &g, const NCConfigVariable &conf); // done
  virtual ~PA_EISMINT_Greenland() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2 &result);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode greenhouse_warming(PetscReal start_year);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2 &result);
protected:
  virtual PetscReal greenhouse_shift(PetscReal t_years, PetscReal dt_years);
  bool do_greenhouse_warming;
  PetscReal greenhouse_warming_start_year;
};

#endif	// __PA_EISMINT_Greenland
