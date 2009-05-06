#!/usr/bin/env python

from netCDF3 import *

input = "../lib/pism_config.nc"

nc = Dataset(input, 'r')

var = nc.variables['pism_config']

print """
/*!
\page config Configuration flags and parameters

\par Notes
Valid configuration flag values are "yes", "true", "on" and "no",
"false", "off" (lowercase only). Note that they have to be inclosed in quotes
in pism_config.cdl.
"""

print """
\section flags Boolean flags

<table style="width: 100%">
<tr>
<td> Flag name </td>
<td> Default value </td>
<td> Description </td>
</tr>
"""

for attr in var.ncattrs():
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "missing")
    except:
      docstring = "missing"

    if type(value) != str:
        continue

    print "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (attr, value, docstring)

print "</table>"

print """
\section params Scalar parameters

<table style="width: 100%">
<tr>
<td> Parameter name </td>
<td> Default value </td>
<td> Description </td>
</tr>
"""

for attr in var.ncattrs():
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "missing")
    except:
      docstring = "missing"

    if type(value) == str:
        continue

    print "<tr><td>%s</td><td>%f</td><td>%s</td></tr>" % (attr, value, docstring)

print """
</table>
*/
"""
