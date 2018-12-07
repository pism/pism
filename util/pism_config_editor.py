#!/usr/bin/env python
# import modules:

# @package pism_config_editor
##
# A script simplifying creating configuration files to use with PISM's -config_override option.
##
# Does not take any command-line options; the only argument is a name of a
# NetCDF configuration file to start from. Run
# \verbatim
# pism_config_editor.py lib/pism_config.nc
# \endverbatim
# to edit lib/pism_config.nc or create a file based on lib/pism_config.nc
##
# \verbatim
# macbook:pism> pism_config_editor.py lib/pism_config.nc
# PISM config file editor: using attributes from 'pism_config' in 'lib/pism_config.nc'.
##
# Please enter a parameter name or hit Return to save your changes.
# You can also hit 'tab' for completions.
# >
# \endverbatim
# Next, start typing the name of a flag or parameter you want to change; hit [tab] to complete:
# \verbatim
# > sta[tab][tab]
# standard_gravity  start_year
# > sta
# \endverbatim
# typing "n[tab][return]" produces:
# \verbatim
# > standard_gravity
##
# Documentation: m s-2; acceleration due to gravity on Earth geoid
# Current value: standard_gravity = 9.8100000000000005
# New value: standard_gravity =
# \endverbatim
# enter the new value (10, for example), press [return]; you would see
# \verbatim
# New value set: standard_gravity = 10.0
##
# List of changes so far:
# standard_gravity = 10.0
##
# Please enter a parameter name or hit Return to save your changes.
# You can also hit 'tab' for completions.
# >
# \endverbatim
##
# Now you can select a different parameter or hit [return] to save to a file:
# \verbatim
# Please enter the file name to save to or hit Return to save to the original file (lib/pism_config.nc).
# > g_equals_10.nc
# \endverbatim
# Next, press [return] if you edited a PISM config file containing \b all the
# parameters or type "pism_overrides[return]" to create a config to use with -config_override.
# \verbatim
# > pism_overrides
# Created variable pism_overrides in g_equals_10.nc.
# Done.
# \endverbatim

import sys
from numpy import double
try:
    import readline
except:
    print("GNU readline library is not available.")
    sys.exit(0)

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


def list_completer(text, state, list):
    """Completes strings from the list 'list'. Skips documenting strings."""
    matches = [x for x in list if (x.startswith(text) and not x.endswith("_doc"))]
    if (state >= len(matches)):
        return None
    else:
        return matches[state]


def edit_attr(dict, attr):
    """Edits an attribute in the dictionary dict."""
    completer = readline.get_completer()
    readline.set_completer(None)
    current_value = dict[attr]
    try:
        print("\n# Documentation: %s" % dict[attr + "_doc"])
    except:
        pass

    print("# Current value: %s = %s" % (attr, str(current_value)))

    while True:
        new_value = eval(input("#     New value: %s = " % attr))

        if new_value == "":
            new_value = current_value
            print("# Using the current value (%s)" % str(current_value))
            break

        try:
            new_value = double(new_value)  # try interpreting as a number
        except:
            pass                # leave as a string

        break

    readline.set_completer(completer)
    return new_value


def main_loop(dict):
    changes = {}
    while True:
        print("\n# Please enter a parameter name or hit Return to save your changes.\n# You can also hit 'tab' for completions.")
        attr = eval(input("> "))
        if attr == "":
            break

        try:
            old_value = dict[attr]
            new_value = edit_attr(dict, attr)
            changes[attr] = new_value
            if (old_value != new_value):
                print("# New value set: %s = %s\n" % (attr, str(new_value)))
        except:
            print("ERROR: attribute '%s' was not found." % attr)

        print("## List of changes so far:")
        for each in list(changes.keys()):
            print("## %s = %s" % (each, str(changes[each])))
    return changes


def read(filename):
    """Reads attributes from a file."""
    try:
        nc = NC(filename)
    except:
        print("ERROR: can't open %s" % filename)
        sys.exit(0)

    names = ['pism_config', 'pism_overrides']
    varname = None
    var = None
    for name in names:
        try:
            var = nc.variables[name]
            varname = name
        except:
            pass

    if var == None:
        print("ERROR: can't find 'pism_config' or 'pism_overrides' in '%s'." % filename)
        sys.exit(0)

    attrs = var.ncattrs()
    dict = {}
    for each in attrs:
        dict[each] = getattr(var, each)
    nc.close()

    return (varname, dict)


def save(dict, changes, default_filename, default_varname):
    """Saves attributes stored in the dictionary changes, adding doc-strings from dict."""
    readline.set_completer(None)

    print("\nPlease enter the file name to save to or hit Return to save to the original file (%s)." % default_filename)
    filename = eval(input("> "))
    if filename == "":
        filename = default_filename

    def varname_completer(text, state):
        names = ['pism_config', 'pism_overrides']
        matches = [x for x in names if x.startswith(text)]

        if state < 2:
            return matches[state]
        else:
            return None

    readline.set_completer(varname_completer)
    print("# Please enter the variable name to use or hit Return to use '%s'." % default_varname)
    varname = eval(input("> "))
    if varname == "":
        varname = default_varname

    try:
        nc = NC(filename, 'a')      # append
    except:
        try:
            nc = NC(filename, 'w', format='NETCDF3_CLASSIC')  # if not found, then create
        except:
            print("ERROR: can't open '%s'." % filename)
            return False

    try:
        var = nc.variables[varname]
    except:
        var = nc.createVariable(varname, 'b')
        print("# Created variable %s in %s." % (varname, filename))

    for each in list(changes.keys()):
        try:
            doc = each + "_doc"
            setattr(var, doc, dict[doc])
        except:
            pass
        setattr(var, each, changes[each])

    nc.close()
    return True


from optparse import OptionParser

parser = OptionParser()

parser.usage = """Run "%prog config.nc"
  to edit config.nc or create a new configuration file (such as an "overrides" file)
  based on config.nc"""
parser.description = "This scrips simplifies creating a customized PISM configuration file."

(options, args) = parser.parse_args()

if (len(args) != 1):
    print("Please specify an input file. Exiting...")
    sys.exit(1)

# Get the input filename:
try:
    filename = args[0]
except:
    sys.exit(0)

# Read attributes:
varname, dict = read(filename)

print("PISM config file editor: using attributes from '%s' in '%s'." % (varname, filename))

# Set up tab completion:


def complete(text, state):
    return list_completer(text, state, list(dict.keys()))


readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

# Process user input:
changes = main_loop(dict)
if changes == {}:
    sys.exit(0)

# Save to a file:
while True:
    result = save(dict, changes, filename, varname)

    if result == True:
        print("Done.")
        break

    print("Do you want to try a different file name? [y/n]")
    answer = eval(input())
    if answer not in ["y", "Y", "yes", "Yes", "YES"]:
        break
