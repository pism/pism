#!/usr/bin/env python
# import modules:
import sys
from numpy import double
try:
    import readline
except:
    print "GNU readline library is not available."
    sys.exit(0)

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

def list_completer(text, state, list):
    """Completes strings from the list 'list'. Skips documenting strings."""
    matches = filter(lambda(x): (x.startswith(text) and not x.endswith("_doc")), list)
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
        print "\n# Documentation: %s" % dict[attr + "_doc"]
    except:
        pass

    print "# Current value: %s = %s" % (attr, str(current_value))

    while True:
        new_value = raw_input("#     New value: %s = " % attr)

        if new_value == "":
            new_value = current_value
            print "# Using the current value (%s)" % str(current_value)
            break

        try:
            new_value = double(new_value) # try interpreting as a number
        except:
            pass                # leave as a string

        break

    readline.set_completer(completer)
    return new_value

def main_loop(dict):
    changes = {}
    while True:
        print "\n# Please enter a parameter name or hit Return to save your changes.\n# You can also hit 'tab' for completions."
        attr = raw_input("> ")
        if attr == "":
            break

        try:
            old_value = dict[attr]
            new_value = edit_attr(dict, attr)
            changes[attr] = new_value
            if (old_value != new_value):
                print "# New value set: %s = %s\n" % (attr, str(new_value))
        except:
            print "ERROR: attribute '%s' was not found." % attr

        print "## List of changes so far:"
        for each in changes.keys():
            print "## %s = %s" % (each, str(changes[each]))
    return changes

def read(filename):
    """Reads attributes from a file."""
    try:
        nc = NC(filename)
    except:
        print "ERROR: can't open %s" % filename
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
        print "ERROR: can't find 'pism_config' or 'pism_overrides' in '%s'." % filename
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

    print "\nPlease enter the file name to save to or hit Return to save to the original file (%s)." % default_filename
    filename = raw_input("> ")
    if filename == "":
        filename = default_filename

    def varname_completer(text, state):
        names = ['pism_config', 'pism_overrides']
        matches = filter(lambda(x): x.startswith(text), names)

        if state < 2:
            return matches[state]
        else:
            return None

    readline.set_completer(varname_completer)
    print "# Please enter the variable name to use or hit Return to use '%s'." % default_varname
    varname = raw_input("> ")
    if varname == "":
        varname = default_varname

    try:
        nc = NC(filename, 'a')      # append
    except:
        try:
            nc = NC(filename, 'w') # if not found, then create
        except:
            print "ERROR: can't open '%s'." % filename
            return False

    try:
        var = nc.variables[varname]
    except:
        var = nc.createVariable(varname, 'b')
        print "# Created variable %s in %s." % (varname, filename)

    for each in changes.keys():
        try:
            doc = each + "_doc"
            setattr(var, doc, dict[doc])
        except:
            pass
        setattr(var, each, changes[each])

    nc.close()
    return True


# Get the input filename:
try:
    filename = sys.argv[1]
except:
    print "Usage: pism_config_editor.py <filename>"
    sys.exit(0)

# Read attributes:
varname, dict = read(filename)

print "PISM config file editor: using attributes from '%s' in '%s'." % (varname, filename)

# Set up tab completion:
def complete(text, state):
    return list_completer(text, state, dict.keys())

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
        print "Done."
        break

    print "Do you want to try a different file name? [y/n]"
    answer = raw_input()
    if answer not in ["y", "Y", "yes", "Yes", "YES"]:
        break
