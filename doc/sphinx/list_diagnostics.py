#!/usr/bin/env python3

import argparse
import json

# file header
header = """.. -*- mode: rst -*-

.. DO NOT EDIT. This file was generated using list_diagnostics.py.
"""

# section start
section = """
.. _{label}:

{title}
{underline}"""

# start of an entry corresponding to a diagnostic
entry_start = """
#. ``{name}``"""

# template used when a diagnostic corresponds to one NetCDF variable
template_single = """
   :Units: {units}
   :Description: {long_name}"""

# template used when a diagnostic corresponds to multiple NetCDF variables
template_many = """
   - ``{var_name}``

     :Units: {units}
     :Description: {long_name}"""

# standard_name line
std_name = """{padding}:Standard name: ``{standard_name}``"""

# comment line
comment = """{padding}:Comment: {comment}"""


def print_diagnostics(title, label, diagnostics):
    print(section.format(label=label, title=title, underline="-" * len(title)))

    for name in sorted(diagnostics.keys()):
        print(entry_start.format(name=name))

        if len(diagnostics[name]) == 1:
            template = template_single
            padding = " " * 3
        else:
            template = template_many
            padding = " " * 5

        for data in diagnostics[name]:
            var_name, units, long_name, standard_name, comment_string = data
            if len(units) == 0:
                units = "---"

            print(template.format(var_name=var_name, units=units, long_name=long_name))

            if len(standard_name) > 0:
                print(std_name.format(padding=padding, standard_name=standard_name))

            if len(comment_string) > 0:
                print(comment.format(padding=padding, comment=comment_string))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = """Generate an RST list documenting PISM's diagnostics."""
    parser.add_argument("FILE", nargs="*")
    options = parser.parse_args()

    spatial = {}
    scalar = {}
    for f in options.FILE:
        with open(f) as f:
            data = json.load(f)
            spatial.update(data["spatial"])
            scalar.update(data["scalar"])

    print(header)

    print_diagnostics("Spatially-variable fields", "sec-extra_vars", spatial)
    print_diagnostics("Scalar time-series", "sec-ts_vars", scalar)
