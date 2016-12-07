/* Classes related to column-wise systems and interpolation in the
 * column (from the storage grid to the fine computational grid).
 *
 * We use these Python wrappers to test C++ code.
 */

%{
#include "verif/tests/exactTestsFG.hh"
#include "verif/tests/exactTestH.h"
#include "verif/tests/exactTestsIJ.h"
#include "verif/tests/exactTestK.h"
#include "verif/tests/exactTestL.hh"
#include "verif/tests/exactTestO.h"
#include "verif/tests/exactTestM.h"
#include "verif/tests/exactTestN.h"
#include "verif/tests/exactTestP.hh"
%}

%include "verif/tests/exactTestsFG.hh"
%include "verif/tests/exactTestH.h"
%include "verif/tests/exactTestsIJ.h"
%include "verif/tests/exactTestK.h"
%include "verif/tests/exactTestL.hh"
%include "verif/tests/exactTestO.h"
%include "verif/tests/exactTestM.h"
%include "verif/tests/exactTestN.h"
%include "verif/tests/exactTestP.hh"
