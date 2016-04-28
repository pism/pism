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

#### Notes
- These flags and parameters are from pism_config.cdl in the src/ directory.
- pism_config.cdl is converted to lib/pism_config.nc in the build process by ncgen.
- An alternate config file foo.nc can be specified at runtime by option "-config foo.nc".
- Values are asked-for by name, by the PISM executable, *when needed*, so if there
  is no request for it then a flag or parameter could be missing and things could still run.
  But the .nc config file must contain any requested values.  The PISM executable
  will terminate if there is no flag or parameter with the requested name.
- Valid boolean flag values are "yes", "true", "on" for TRUE and "no",
  "false", "off" for FALSE.  Lowercase only.  They have to be enclosed in quotes
  in pism_config.cdl.

#### To create and use an alternate config file:

#### Method 1, by using util/pism_config_editor.py:
- Explained at pism_config_editor.  It may be the easiest way!

#### Method 2, by editing a .cdl text file:
- Make a copy of src/pism_config.cdl:
~~~
    cp src/pism_config.cdl myconfig.cdl
~~~
- Edit the text file myconfig.cdl to have the values you want.  Generally values of parameters
  can be changed but it is dangerous to remove them entirely.  If you are building a derived class
  of IceModel then you might add new values.
- Create a new configuration .nc file with ncgen, and do your run with the new values.  For example:
~~~
    ncgen -o myconfig.nc myconfig.cdl
    pismr -config myconfig.nc -i mydata.nc -bootstrap -Mx 101 -My 101 -Mz 101 -Lz 4000 -y 100 -o start.nc
    pismr -config myconfig.nc -i start.nc -y 10000 -o end.nc
~~~
  (Runtime option "-verbose 4" will report back your values as the PISM
  executable starts.)

#### Method 3, using a netCDF Operator (NCO):
- This illustration changes the Clausius-Clapeyron constant from its default value to
  9.7008e-8 K Pa-1.  First you make a copy of lib/pism_config.nc, assuming the
  PISM source is built.  Then make your modification of the desired attribute
  (of the only variable in myconfig.nc, namely `pism_config`), using ncatted
  (see <a href="http://nco.sourceforge.net/">NCO homepage</a>).  Then view your handiwork
  with ncdump:
~~~
    cp lib/pism_config.nc myconfig.nc
    ncatted -a ice.beta_Clausius_Clapeyron,pism_config,m,d,9.7008e-8 myconfig.nc
    ncdump -h myconfig.nc | grep ice.beta_Clausius_Clapeyron
~~~
  Now run with the new values as before:
~~~
    pismr -config myconfig.nc -i mydata.nc -bootstrap -Mx 101 -My 101 -Mz 101 -Lz 4000 -y 100 -o start.nc
    ...
~~~
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
    print '<table style="width: 100%">\n'

    print_row([titles[0]] + ["<b>{}</b>".format(t) for t in titles[1:]])

def print_row(data):
    "Print a row of a table."
    print '<tr>'
    print '<td class="indexkey">%s</td>' % data[0]
    for d in data[1:]:
        print '<td class="indexvalue">%s</td>' % d
    print '</tr>'

def print_footer():
    print "</table>"

def get_units(var, name):
    "Get the units string for a parameter."
    return getattr(var, name + "_units", "[none]")

def get_docstring(var, name):
    "Get the documenting string for a parameter."
    return getattr(var, name + "_doc", "[missing]")

def is_special(var, name):
    "Check if the name is 'special' and should not be included."

    if name == "long_name":
        return True

    for n in ["_doc", "_units", "_type", "_option", "_choices"]:
        if name.endswith(n):
            return True

    return False

def is_boolean(var, name):
    "Check if a name corresponds to a boolean flag."
    return (not is_special(var, name)) and (getattr(var, name) in allowed_boolean_values)

def is_number(var, name):
    "Check if a name corresponds to a scalar parameter."
    return (not is_special(var, name)) and isinstance(getattr(var, name), (int, float))

def is_string(var, name):
    "Check if a name corresponds to a string."
    return (not is_special(var, name)) and isinstance(getattr(var, name), (str, unicode)) and (not is_boolean(var, name))

def print_parameters(var, parameter_list, transform=lambda x: "\"%s\"" % x, use_units=False, use_choices=False):
    "Print table rows corresponding to parameters in a list."
    for name in sorted(parameter_list):
        option = "`-" + getattr(var, name + "_option", name) + "`"
        if use_units:
            print_row((name,
                       transform(getattr(var, name)),
                       get_units(var, name),
                       option,
                       get_docstring(var, name)))
        elif use_choices:
            # replace comma with comma-space to allow wrapping
            choices = "`" + getattr(var, name + "_choices", "---").replace(",", ", ") + "`"
            print_row((name,
                       transform(getattr(var, name)),
                       choices,
                       option,
                       get_docstring(var, name)))
        else:
            print_row((name,
                       transform(getattr(var, name)),
                       option,
                       get_docstring(var, name)))

def number_to_string(number):
    "Format a number as a string. Use scientific notation for large and small numbers, remove '.0' from integers."
    if abs(number) >= 1e7:
        return '%e' % number # use scientific notation if a number is big
    elif int(number) == number:
        return '%d' % int(number) # remove zeros after the decimal point
    elif abs(number) <= 1e-5:
        return '%e' % number # use scientific notation if small (and not zero; previous case)
    else:
        return '%f' % number

def print_booleans(var):
    print "@section flags Boolean flags"
    print_header(("Flag name", "Value", "Command-line option", "Description"))
    print_parameters(var,
                     filter(lambda name: is_boolean(var, name),
                            var.ncattrs()))
    print_footer()

def print_scalars(var):
    print "@section params Scalar parameters"
    print_header(("Parameter name", "Value", "Units", "Command-line option", "Description"))
    print_parameters(var,
                     filter(lambda name: is_number(var, name),
                            var.ncattrs()),
                     number_to_string,
                     use_units=True)
    print_footer()

def print_strings(var):
    print "@section strings String parameters"
    print_header(("Parameter name", "Value", "Allowed values", "Command-line option", "Description"))
    print_parameters(var,
                     filter(lambda name: is_string(var, name),
                            var.ncattrs()), use_choices=True)
    print_footer()

if __name__ == "__main__":
    print intro

    nc = netCDF4.Dataset(input, 'r')

    var = nc.variables['pism_config']

    print_booleans(var)

    print_scalars(var)

    print_strings(var)
