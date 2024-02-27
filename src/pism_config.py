#!/usr/bin/env python3

import netCDF4
import numpy
import json
import sys
import argparse

suffixes = ["choices", "doc", "option", "type", "units"]

def special(name):
    "Return true if 'name' is a 'special' parameter, false otherwise."
    for suffix in suffixes:
        if name.endswith("_" + suffix) or name == "long_name":
            return True
    return False

def nested_dict(config):
    """Create a nested Python dictionary corresponding to parameters stored in a "flat"
    dictionary 'config'.

    """
    result = {}
    for name in sorted(config.keys()):
        # skip "special" attributes
        if special(name):
            continue

        # create the nested structure
        leaf = result
        for t in name.split("."):
            if not t in leaf.keys():
                leaf[t] = {}
            leaf = leaf[t]

        # Convert numpy.int32 to int (json.dump() cannot handle np.int32)
        attr = config[name]
        if isinstance(attr, numpy.int32):
            attr = int(attr)

        # convert flags to Booleans
        if config[name + "_type"] == "flag":
            attr = True if attr == "yes" else False

        # store the value
        leaf["value"] = attr

        # store metadata
        for s in suffixes:
            try:
                leaf[s] = config[name + "_" + s]
            except:
                pass

    return result

def flat_dict(config, path=None):
    """Convert a nested Python dictionary 'config' describing PISM's configuration parameters
    into a flat one.

    """
    result = {}
    for key in sorted(config.keys()):
        if isinstance(config[key], dict):
            result.update(flat_dict(config[key], path + "." + key if path else key))
        else:
            if key == "value":
                result[path] = config[key]
            elif key in suffixes:
                result[path + "_" + key] = config[key]
            else:
                raise ValueError("unknown leaf: {}".format(key))
    return result

def netcdf_to_flat_dict(variable):
    "Create a flat dictionary containing all attributes of a NetCDF variable."
    return {k : getattr(variable, k) for k in variable.ncattrs()}

def cdl(config):
    """Return the NetCDF CDL description of a PISM configuration database stored in a flat
    Python dictionary.
    """

    header = """netcdf pism_config {
    variables:
    byte pism_config;

"""

    footer = """    pism_config:long_name = "PISM configuration flags and parameters.";
    pism_config:long_name_doc = "The 'long_name' attribute is required by CF conventions. It is not used by PISM itself.";
}
"""

    result = header
    for key in sorted(config.keys()):
        if not special(key):
            if config[key + "_type"] in ["string", "keyword", "flag"]:
                result += '    pism_config:{} = "{}";\n'.format(key, config[key])
            else:
                result += '    pism_config:{} = {};\n'.format(key, config[key])
            for s in suffixes:
                try:
                    value = config[key + "_" + s]
                    if s in ["doc", "units"]:
                        value = value.replace("\\", "\\\\") # escape backslashes in formulas
                        value = value.replace('"', '\\"')   # escape double quotes
                    result += '    pism_config:{}_{} = "{}";\n'.format(key, s, value)
                except:
                    pass
            result += "\n"

    result += footer

    return result

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("FILE", nargs=1)
    parser.add_argument("--json", dest="json", action="store_true",
                        help="output JSON", default=False)
    options = parser.parse_args()

    config = netcdf_to_flat_dict(netCDF4.Dataset(options.FILE[0]).variables["pism_config"])

    if options.json:
        print(json.dumps(nested_dict(config), indent=4, sort_keys=True))
    else:
        print(cdl(config))
