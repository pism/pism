#!/usr/bin/env python
# Checks the structure of a PISM configuration file. (This is used as a regression test.)


def is_special(name):
    "Check if the name is 'special' and should not be included."

    if name == "long_name":
        return True

    for n in ["_doc", "_units", "_type", "_option", "_choices"]:
        if name.endswith(n):
            return True

    return False


import netCDF4
import sys
import numpy as np

config = netCDF4.Dataset(sys.argv[1])

pism_config = config.variables['pism_config']

attrs = pism_config.ncattrs()

for a in attrs:
    if is_special(a):
        continue

    attr_value = getattr(pism_config, a)

    if (a + "_doc") not in attrs:
        print("Attribute {} is not documented".format(a))
        sys.exit(1)

    if (a + "_type") not in attrs:
        print("Attribute {} does not have a type".format(a))
        sys.exit(1)

    attr_type = getattr(pism_config, a + "_type")

    if attr_type in ["scalar", "integer"] and (a + "_units") not in attrs:
        print("Attribute {} is a number, but it does not have units".format(a))
        sys.exit(1)

    if attr_type == "boolean" and attr_value not in ["yes", "no", "true", "false", "on", "off"]:
        print("Attribute {} is a boolean, but its value is {}".format(a, attr_value))
        sys.exit(1)

    if attr_type == "keyword" and (a + "_choices") not in attrs:
        print("Attribute {} is a keyword, but {}_choices is not set".format(a, a))
        sys.exit(1)

    if attr_type in ["scalar", "integer"] and not isinstance(attr_value, (int, float, np.number)):
        print("Attribute {} is a number, but its value is not".format(a))
        print(type(attr_value))
        sys.exit(1)

    if attr_type not in ["scalar", "integer"] and isinstance(attr_value, (int, float)):
        print("Attribute {} is not a number, but its value is".format(a))
