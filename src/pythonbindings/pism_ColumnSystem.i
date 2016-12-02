/* Classes related to column-wise systems and interpolation in the
 * column (from the storage grid to the fine computational grid).
 *
 * We use these Python wrappers to test C++ code.
 */

%{
#include "base/util/ColumnInterpolation.hh"
%}

%include "base/util/ColumnInterpolation.hh"
