#!/usr/bin/env python3

"""This script prints the 'nested' JSON representation of a PISM configuration contained
in a provided file (in the 'pism_config' variable). It is not used anywhere at the moment,
but may be useful in the future."""

import netCDF4
import json
import argparse

suffixes = ["choices", "doc", "option", "type", "units", "valid_min", "valid_max"]

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

        attr = config[name]

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

def int32_converter(obj):
    "Convert int32 to int to make them JSON-serializable"
    try:
        return int(obj)
    except:
        raise TypeError(f"Cannot convert {obj} to int (type: {obj.__class__.__name__})")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("FILE", nargs=1)
    options = parser.parse_args()

    pism_config = netCDF4.Dataset(options.FILE[0]).variables["pism_config"]
    attributes = {k : getattr(pism_config, k) for k in pism_config.ncattrs()}
    config = nested_dict(attributes)

    json_config = json.dumps(config, indent=2, sort_keys=True, default=int32_converter, allow_nan=False)

    print(json_config)
