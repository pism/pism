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

std_name_template = """{padding}:Standard name: ``{standard_name}``"""

comment_template = """{padding}:Comment: {comment}"""

template_many = """
   - ``{var_name}``

     :Units: {units}
     :Description: {long_name}"""

header = """
.. DO NOT EDIT. This file was generated using list_diagnostics.py.

   To update this list, run PISM with the option -list_diagnostics_json, making sure that
   you enabled all the sub-models that provide diagnostics (scalar and
   spatially-variable).

   At the time of writing this includes

   - Evolving geometry (without -no_mass)
   - The bedrock thermal layer model (-Mbz X for X > 1)
   - An energy-conservation model (-energy enthalpy or -energy cold)
   - A bed deformation model (-bed_def iso or -bed_def lc)
   - The PDD surface model (-surface pdd)
   - An atmosphere model using a cosine yearly cycle (e.g. -atmosphere searise_greenland)
   - The hybric stress balance model (-stress_balance ssa+sia -ssa_method fd)
   - The Mohr-Coulomb basal yield stress model (-yield_stress mohr_coulomb)
   - The "routing" subglacial hydrology model (-hydrology routing)
   - All supported calving mechanisms (-calving eigen_calving,thickness_calving,frontal_melt,ocean_kill,vonmises_calving)

   The "distributed" hydrology model does add more diagnostics, but it is not supported at
   this point (not documented in the User's Manual).

   See the Makefile in this directory.

.. _sec-diagnostics-list:

List of PISM's diagnostics
==========================

.. contents::

.. include:: diagnostics-comments.txt

"""

def print_diagnostics(diagnostics):

    def print_some(title, label, diagnostics):

        print_heading(title, label, "-")

        for name in sorted(diagnostics.keys()):
            print(entry_start.format(name=name))

            if len(diagnostics[name]) == 1:
                template = template_single
                padding = " " * 3
            else:
                template = template_many
                padding = " " * 5

            for data in diagnostics[name]:
                var_name, units, long_name, standard_name, comment = data
                if len(units) == 0:
                    units = "---"

                print(template.format(var_name=var_name,
                                      units=units,
                                      long_name=long_name))

                if len(standard_name) > 0:
                    print(std_name_template.format(padding=padding,
                                                   standard_name=standard_name))

                if len(comment) > 0:
                    print(comment_template.format(padding=padding,
                                                  comment=comment))

    print(header)

    print_some("Spatially-variable fields", "sec-extra_vars", diagnostics["spatial"])
    print_some("Scalar time-series",  "sec-ts_vars", diagnostics["scalar"])

with open(options.FILE[0]) as f:
    print_diagnostics(json.load(f))
