/* Classes related to column-wise systems and interpolation in the
 * column (from the storage grid to the fine computational grid).
 *
 * We use these Python wrappers to test C++ code.
 */

%{
#include "enthSystem.hh"
#include "ColumnInterpolation.hh"
%}

/* wrap the enthalpy solver to make testing easier */
%include "columnSystem.hh"
%rename(get_lambda) pism::enthSystemCtx::lambda;
%include "enthSystem.hh"

%include "ColumnInterpolation.hh"
