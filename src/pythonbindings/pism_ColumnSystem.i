/* Classes related to column-wise systems and interpolation in the
 * column (from the storage grid to the fine computational grid).
 *
 * We use these Python wrappers to test C++ code.
 */

%{
#include "util/ColumnInterpolation.hh"
%}

%include "util/ColumnInterpolation.hh"
