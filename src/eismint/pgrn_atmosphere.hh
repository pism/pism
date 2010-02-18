#ifndef __PA_EISMINT_Greenland
#define __PA_EISMINT_Greenland

#include "../coupler/PISMAtmosphere.hh"

class PA_EISMINT_Greenland : public PA_Parameterized_Temperature {
public:
  PA_EISMINT_Greenland(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PA_EISMINT_Greenland() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
protected:
  virtual PetscReal greenhouse_shift(PetscReal t_years, PetscReal dt_years);
  bool do_greenhouse_warming;
  PetscReal greenhouse_warming_start_year;
  IceModelVec2 *lat, *surfelev;
};

#endif	// __PA_EISMINT_Greenland
