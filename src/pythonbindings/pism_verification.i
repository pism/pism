/* Classes related to column-wise systems and interpolation in the
 * column (from the storage grid to the fine computational grid).
 *
 * We use these Python wrappers to test C++ code.
 */

%{
#include "verification/tests/exactTestsABCD.h"
#include "verification/tests/exactTestsFG.hh"
#include "verification/tests/exactTestH.h"
#include "verification/tests/exactTestsIJ.h"
#include "verification/tests/exactTestK.h"
#include "verification/tests/exactTestL.hh"
#include "verification/tests/exactTestO.h"
#include "verification/tests/exactTestM.h"
#include "verification/tests/exactTestN.h"
#include "verification/tests/exactTestP.hh"
%}

%include "verification/tests/exactTestsABCD.h"
%include "verification/tests/exactTestsFG.hh"
%include "verification/tests/exactTestH.h"
%include "verification/tests/exactTestsIJ.h"
%include "verification/tests/exactTestK.h"
%include "verification/tests/exactTestL.hh"
%include "verification/tests/exactTestO.h"
%include "verification/tests/exactTestM.h"
%include "verification/tests/exactTestN.h"
%include "verification/tests/exactTestP.hh"
