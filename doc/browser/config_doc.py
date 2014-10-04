#!/usr/bin/env python

intro = """
Configuration flags and parameters {#config}
=============

[TOC]

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

try:
    import netCDF4
except:
    print "ERROR: netCDF4 is required!"
    import sys
    sys.exit(1)

input = "pism_config.nc"

allowed_boolean_values = ["yes", "true", "on", "no", "false", "off"]

def print_header(titles):
    "Print the header of a table. 'title' is a tuple with three strings, the titles for table columns."
    print '<table style="width: 100%">\n' + \
    '<tr> <td class="indexkey"> %s </td> <td class="indexvalue"> <b> %s </b> </td> <td class="indexvalue"> <b> %s </b> </td> </tr>' % titles

def print_row(name, value, docstring):
    "Print a row of a table."
    print '<tr><td class="indexkey">%s</td><td class="indexvalue"> %s </td><td class="indexvalue">%s</td></tr>' % (name, value, docstring)

def print_footer():
    print "</table>"

def get_docstring(name):
    "Get the documenting string for a parameter."
    global var
    return getattr(var, name + "_doc", "[missing]")

def is_special(name):
    "Check if the name is 'special' and should not be included."
    global var
    return (name in ("long_name")) or name.endswith("_doc")

def is_boolean(name):
    "Check if a name corresponds to a boolean flag."
    global var
    return (not is_special(name)) and (getattr(var, name) in allowed_boolean_values)

def is_number(name):
    "Check if a name corresponds to a scalar parameter."
    global var
    return (not is_special(name)) and isinstance(getattr(var, name), (int, float))

def is_string(name):
    "Check if a name corresponds to a string."
    global var
    return (not is_special(name)) and isinstance(getattr(var, name), (str, unicode)) and (not is_boolean(name))

def print_parameters(parameter_list, transform=lambda x: "\"%s\"" % x):
    "Print table rows corresponding to parameters in a list."
    global var
    for name in sorted(parameter_list):
        print_row(name, transform(getattr(var, name)), get_docstring(name))

def number_to_string(number):
    "Format a number as a string. Use scientific notation for large and small numbers, remove '.0' from integers."
    if (abs(number) >= 1e7):
        return '%e' % number # use scientific notation if a number is big
    elif (int(number) == number):
        return '%d' % int(number) # remove zeros after the decimal point
    elif (abs(number) <= 1e-5):
        return '%e' % number # use scientific notation if small (and not zero; previous case)
    else:
        return '%f' % number

def print_booleans():
    global var
    print "\section flags Boolean flags"
    print_header(("Flag name", "Default value", "Description"))
    print_parameters(filter(is_boolean, var.ncattrs()))
    print_footer()

def print_scalars():
    global var
    print "\section params Scalar parameters"
    print_header(("Parameter name", "Default value", "Description"))
    print_parameters(filter(is_number, var.ncattrs()), number_to_string)
    print_footer()

def print_strings():
    global var
    print "\section strings String parameters"
    print_header(("Parameter name", "Default value", "Description"))
    print_parameters(filter(is_string, var.ncattrs()))
    print_footer()

if __name__ == "__main__":
    print intro

    nc = netCDF4.Dataset(input, 'r')

    var = nc.variables['pism_config']

    print_booleans()

    print_scalars()

    print_strings()
