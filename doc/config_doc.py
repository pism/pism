#!/usr/bin/env python

from netCDF3 import *

input = "../lib/pism_config.nc"

nc = Dataset(input, 'r')

var = nc.variables['pism_config']

print """
/*!
\page config Configuration flags and parameters

\par Notes
- These flags and parameters are from pism_config.cdl in the src/ directory.
- pism_config.cdl is converted to lib/pism_config.nc in the build process by ncgen.
- An alternate config file foo.nc can be specified at runtime by option "-config foo.nc".
- Values are asked-for by name, by the PISM executable, \e when \e needed, so if there
  is no request for it then a flag or parameter could be missing and things could still run.
- But the .nc config file must contain any requested values.  The PISM executable
  will terminate if there is no flag or parameter with the requested name.
- Valid boolean flag values are "yes", "true", "on" for TRUE and "no",
  "false", "off" for FALSE.  Lowercase only.  They have to be enclosed in quotes
  in pism_config.cdl.
  
\par One way to create and use an alternate config file:
- Make a copy of src/pism_config.cdl:
  \code
    cp src/pism_config.cdl myconfig.cdl
  \endcode
- Edit the text file myconfig.cdl to have the values you want.  Generally values of parameters
  can be changed but it is dangerous to remove them entirely.  If you are building a derived class
  of IceModel then you might add new values.
- Create a new configuration .nc file with ncgen:
  \code
    ncgen -o myconfig.nc myconfig.cdl
  \endcode
- Do your run with the new values, for example:
  \code
    pismr -config myconfig.nc -boot_from mydata.nc -Mx 101 -My 101 -Mz 101 -Lz 4000 -y 100 -o start.nc
    pismr -config myconfig.nc -i start.nc -y 10000 -o end.nc
  \endcode
- Runtime option "-verbose 4" will report back your values as the PISM
  executable starts.
"""

print """
\section flags Boolean flags
<table style="width: 100%">
<tr> <td> <b> Flag name </b> </td> <td> <b> Default value </b> </td> <td> <b> Description </b> </td> </tr>"""

for attr in var.ncattrs():
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "[missing]")
    except:
      docstring = "[missing]"

    if type(value) != str:
        continue

    print "<tr><td>%s</td><td>\"%s\"</td><td>%s</td></tr>" % (attr, value, docstring)

print "</table>"

print """
\section params Scalar parameters
<table style="width: 100%">
<tr> <td> <b> Flag name </b> </td> <td> <b> Default value </b> </td> <td> <b> Description </b> </td> </tr>"""

for attr in var.ncattrs():
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "[missing]")
    except:
      docstring = "[missing]"

    if type(value) == str:
        continue

    print "<tr><td>%s</td>" % attr
    
    if (abs(value) >= 1e7):
        print "<td>%e</td>" % value, # use scientific notation if a number is big
    elif (int(value) == value):
        print "<td>%d</td>" % int(value), # remove zeros after the decimal point
    elif (abs(value) <= 1e-5):
        print "<td>%e</td>" % value, # use scientific notation if small (and not zero; prev case)
    else:
        print "<td>%f</td>" % value,
    
    print "<td>%s</td></tr>" % docstring

print "</table> */"
