#!/usr/bin/env python

import json, argparse

parser = argparse.ArgumentParser()
parser.description = '''Generate an RST list documenting PISM's diagnostics.'''
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()

def print_heading(title, label, decoration):
    heading = """
.. _{label}:

{title}
{underline}"""
    print(heading.format(label=label, title=title, underline=decoration * len(title)))

entry_start = """
#. ``{name}``"""

template_single = """
   :Units: {units}
   :Description: {long_name}"""

std_name_single = """   :Standard name: ``{standard_name}``"""

comment_single = """   :Comment: {comment}"""

template_many = """
   - ``{var_name}``

     :Units: {units}
     :Description: {long_name}"""

comment_many = """     :Comment: {comment}"""

std_name_many = """     :Standard name: ``{standard_name}``"""


header = """
.. DO NOT EDIT. This file was generated using list_diagnostics.py.

.. contents::"""

def print_diagnostics(diagnostics):

    def print_some(title, label, diagnostics):

        print_heading(title, label, "-")

        for name in sorted(diagnostics.keys()):
            print(entry_start.format(name=name))

            if len(diagnostics[name]) == 1:
                template = template_single
                std_name_template = std_name_single
                comment_template = comment_single
            else:
                template = template_many
                std_name_template = std_name_many
                comment_template = comment_many

            for data in diagnostics[name]:
                var_name, units, long_name, standard_name, comment = data
                if len(units) == 0:
                    units = "---"

                print(template.format(var_name=var_name,
                                      units=units,
                                      long_name=long_name))

                if len(standard_name) > 0:
                    print(std_name_template.format(standard_name=standard_name))

                if len(comment) > 0:
                    print(comment_template.format(comment=comment))

    print_heading("List of PISM's diagnostics", "sec-diagnostics-list", "=")
    print(header)

    print_some("Spatially-variable fields", "sec-extra_vars", diagnostics["spatial"])
    print_some("Scalar time-series",  "sec-ts_vars", diagnostics["scalar"])

with open(options.FILE[0]) as f:
    print_diagnostics(json.load(f))
