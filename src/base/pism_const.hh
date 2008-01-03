#ifndef __pism_const_hh
#define __pism_const_hh

const PetscScalar gasConst_R = 8.31441;      // J/(mol K)    Gas Constant
const PetscScalar grav       = 9.81;         // m/s^2        acceleration of gravity
const PetscScalar secpera    = 3.1556926e7;
const PetscScalar pi         = 3.14159265358979;

const PetscInt HISTORY_STRING_LENGTH = 32768; // 32KiB ought to be enough

#endif
