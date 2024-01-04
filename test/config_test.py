#!/usr/bin/env python3
# Checks the structure of a PISM configuration file. (This is used as a regression test.)


def is_special(name):
    "Check if the name is 'special' and should not be included."

    if name == "long_name":
        return True

    for n in ["_doc", "_units", "_type", "_option", "_choices"]:
        if name.endswith(n):
            return True

    return False


import xarray as xr
import sys
import numpy as np

config = xr.open_dataset(sys.argv[1])

pism_config = config.variables['pism_config']

attrs = pism_config.attrs

for a in attrs:
    if is_special(a):
        continue

    attr_value = pism_config.attrs[a]

    if (a + "_doc") not in attrs:
        print(f"Attribute {a} is not documented")
        sys.exit(1)

    if (a + "_type") not in attrs:
        print(f"Attribute {a} does not have a type")
        sys.exit(1)

    attr_type = pism_config.attrs[a + "_type"]

    if attr_type in ["number", "integer"] and (a + "_units") not in attrs:
        print("Attribute {} is a number, but it does not have units".format(a))
        sys.exit(1)

    if attr_type == "flag" and attr_value not in ["yes", "no", "true", "false", "on", "off"]:
        print("Attribute {} is a flag, but its value is {}".format(a, attr_value))
        sys.exit(1)

    if attr_type == "keyword" and (a + "_choices") not in attrs:
        print("Attribute {} is a keyword, but {}_choices is not set".format(a, a))
        sys.exit(1)

    if attr_type in ["number", "integer"] and not isinstance(attr_value, (int, float, np.number)):
        print("Attribute {} is a number, but its value is not".format(a))
        print(type(attr_value))
        sys.exit(1)

    if attr_type not in ["number", "integer"] and isinstance(attr_value, (int, float)):
        print("Attribute {} is not a number, but its value is".format(a))
