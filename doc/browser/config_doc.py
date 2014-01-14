#!/usr/bin/env python
from numpy import sort
# try different netCDF modules
try:
    from netCDF4 import *
except:
    try:
        from netCDF3 import *
    except:
        print "ERROR: Neither netCDF3 nor netCDF4 is installed!"
        import sys
        sys.exit(1)

input = "pism_config.nc"

nc = Dataset(input, 'r')

var = nc.variables['pism_config']

print """
# Configuration flags and parameters {#config}

\htmlonly
<p style="text-align: right">
With four parameters I can fit an elephant, and with five I can make him wiggle his trunk.<br><br>
John von Neumann
</p>
\endhtmlonly

\par Notes
- These flags and parameters are from pism_config.cdl in the src/ directory.
- pism_config.cdl is converted to lib/pism_config.nc in the build process by ncgen.
- An alternate config file foo.nc can be specified at runtime by option "-config foo.nc".
- Values are asked-for by name, by the PISM executable, \e when \e needed, so if there
  is no request for it then a flag or parameter could be missing and things could still run.
  But the .nc config file must contain any requested values.  The PISM executable
  will terminate if there is no flag or parameter with the requested name.
- Valid boolean flag values are "yes", "true", "on" for TRUE and "no",
  "false", "off" for FALSE.  Lowercase only.  They have to be enclosed in quotes
  in pism_config.cdl.

\par To create and use an alternate config file:

\par Method 1, by using util/pism_config_editor.py:
- Explained at pism_config_editor.  It may be the easiest way!

\par Method 2, by editing a .cdl text file:
- Make a copy of src/pism_config.cdl:
  \code
    cp src/pism_config.cdl myconfig.cdl
  \endcode
- Edit the text file myconfig.cdl to have the values you want.  Generally values of parameters
  can be changed but it is dangerous to remove them entirely.  If you are building a derived class
  of IceModel then you might add new values.
- Create a new configuration .nc file with ncgen, and do your run with the new values.  For example:
  \code
    ncgen -o myconfig.nc myconfig.cdl
    pismr -config myconfig.nc -boot_file mydata.nc -Mx 101 -My 101 -Mz 101 -Lz 4000 -y 100 -o start.nc
    pismr -config myconfig.nc -i start.nc -y 10000 -o end.nc
  \endcode
  (Runtime option "-verbose 4" will report back your values as the PISM
  executable starts.)

\par Method 3, using a netCDF Operator (NCO):
- This illustration changes the Clausius-Clapeyron constant from its default value to
  9.7008e-8 K Pa-1.  First you make a copy of lib/pism_config.nc, assuming the
  PISM source is built.  Then make your modification of the desired attribute
  (of the only variable in myconfig.nc, namely \c pism_config), using ncatted
  (see <a href="http://nco.sourceforge.net/">NCO homepage</a>).  Then view your handiwork
  with ncdump:
  \code
    cp lib/pism_config.nc myconfig.nc
    ncatted -a beta_CC,pism_config,m,d,9.7008e-8 myconfig.nc
    ncdump -h myconfig.nc | grep beta_CC
  \endcode
  Now run with the new values as before:
  \code
    pismr -config myconfig.nc -boot_file mydata.nc -Mx 101 -My 101 -Mz 101 -Lz 4000 -y 100 -o start.nc
    ...
  \endcode
"""

print """
\section flags Boolean flags
<table style="width: 100%">
<tr> <td class="indexkey"> Flag name </td> <td class="indexvalue"> <b> Default value </b> </td> <td class="indexvalue"> <b> Description </b> </td> </tr>"""

for attr in sort(var.ncattrs()):
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "[missing]")
    except:
      docstring = "[missing]"

    # ignore anything that does not represent a boolean:
    if (value not in ["yes", "true", "on", "no", "false", "off"]):
        continue

    print '<tr><td class="indexkey">%s</td><td class="indexvalue">\"%s\"</td><td class="indexvalue">%s</td></tr>' % (attr, value, docstring)

print "</table>"

print """
\section params Scalar parameters
<table style="width: 100%">
<tr> <td class="indexkey"> <b>Parameter name </b> </td> <td class="indexvalue"> <b> Default value </b> </td> <td class="indexvalue"> <b> Description </b> </td> </tr>"""

for attr in sort(var.ncattrs()):
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "[missing]")
    except:
      docstring = "[missing]"

    if isinstance(value, (str, unicode)):
        continue

    print '<tr><td class="indexkey">%s</td>' % attr

    if (abs(value) >= 1e7):
        print '<td class="indexvalue">%e</td>' % value, # use scientific notation if a number is big
    elif (int(value) == value):
        print '<td class="indexvalue">%d</td>' % int(value), # remove zeros after the decimal point
    elif (abs(value) <= 1e-5):
        print '<td class="indexvalue">%e</td>' % value, # use scientific notation if small (and not zero; prev case)
    else:
        print '<td class="indexvalue">%f</td>' % value,

    print unicode('<td class="indexvalue">%s</td></tr>' % docstring).encode('utf8')

print "</table>"

print """
\section strings String parameters
<table style="width: 100%">
<tr> <td class="indexkey"> <b> Parameter name </b> </td> <td class="indexvalue"> <b> Default value </b> </td> <td class="indexvalue"> <b> Description </b> </td> </tr>"""

for attr in sort(var.ncattrs()):
    if attr.endswith("_doc"):
        continue

    value = getattr(var, attr)
    try:
      docstring = getattr(var, attr + "_doc", "[missing]")
    except:
      docstring = "[missing]"

    # ignore non-strings and strings representing booleans
    if (type(value) != str) or (value in ["yes", "true", "on", "no", "false", "off"]):
        continue

    print '<tr><td class="indexkey">%s</td><td class="indexvalue">\'%s\'</td><td class="indexvalue">%s</td></tr>' % (attr, value, docstring)

print "</table>"
