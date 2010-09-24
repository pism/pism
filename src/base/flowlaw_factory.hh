#ifndef __flowlaw_factory
#define __flowlaw_factory

#include "flowlaws.hh"

#define ICE_CUSTOM  "custom"        /* Plain isothermal Glen with customizable parameters */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (PolyThermalGPBLDIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_HYBRID  "hybrid"        /* Goldsby-Kohlstedt for SIA, PB for SSA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

class IceFlowLawFactory {
public:
  IceFlowLawFactory(MPI_Comm,const char prefix[], const NCConfigVariable &conf);
  ~IceFlowLawFactory();
  PetscErrorCode setType(const char[]);
  PetscErrorCode setFromOptions();
  PetscErrorCode registerType(const char[],
		  PetscErrorCode(*)(MPI_Comm,const char[], const NCConfigVariable &,IceFlowLaw **));
  PetscErrorCode create(IceFlowLaw **);
private:
  PetscErrorCode registerAll();
private:
  MPI_Comm comm;
  char prefix[256],type_name[256];
  PetscFList type_list;
  const NCConfigVariable &config;
};


#endif  // __flowlaw_factory
