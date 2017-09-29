#!/usr/bin/env python

import netCDF4, argparse

def is_special(name):
    "Check if the name is 'special' and should not be included."

    if name == "long_name":
        return True

    for n in ["_doc", "_units", "_type", "_option", "_choices"]:
        if name.endswith(n):
            return True

    return False

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

def entry(name, value, option, description, choices=None):

    if len(value) == 0:
        value = "*no default*"

    if option is not None:
        option = ":opt:`{}`".format(option)
    else:
        option = "*no short option*"

    if value.find(",") != -1:
        # A comma-separated list
        value = "``{}``".format(value.replace(",", ", "))

    if choices is not None:
        # This is always a comma-separated list.
        choices = "``{}``".format(choices.replace(",", ", "))

    template = """
#. :config:`{name}`

   :Value: {value}"""

    if choices is not None:
        template += """
   :Choices: {choices}"""

    template += """
   :Option: {option}
   :Description: {description}"""

    return template.format(name=name,
                           value=value,
                           choices=choices,
                           option=option,
                           description=description)

def print_string(var, name):
    print(entry(name,
                value(var, name),
                option(var, name),
                doc(var, name)))

def print_scalar(var, name):
    V = "{} ({})".format(number_to_string(value(var, name)), units(var, name))
    print(entry(name, V, option(var, name), doc(var, name)))

def print_integer(var, name):
    print(entry(name,
                str(value(var, name)),
                option(var, name),
                doc(var, name)))

def print_keyword(var, name):
    print(entry(name,
                value(var, name),
                option(var, name),
                doc(var, name),
                choices=value(var, name + "_choices")))

def value(var, name):
    return var.getncattr(name)

def units(var, name):
    return value(var, name + "_units")

def doc(var, name):
    return value(var, name + "_doc")

def option(var, name):
    try:
        return "-" + value(var, name + "_option")
    except AttributeError:
        return None

printers = {"string" : print_string,
            "scalar" : print_scalar,
            "integer" : print_integer,
            "boolean" : print_string,
            "keyword" : print_keyword}

header = """
.. DO NOT EDIT: This file was automatically generated using config_parameters.py. Edit src/pism_config.cdl instead.

.. include:: ../../global.rst

.. _sec-parameter-list:

List of configuration parameters
================================

Each parameter can be set using the command-line option consisting of a dash followed by
the parameter name. For example,

.. code-block:: none

   -constants.standard_gravity 10

sets the acceleration due to gravity (parameter :config:`constants.standard_gravity`) to
`10`.

"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("FILE", nargs=1)
    options = parser.parse_args()

    f = netCDF4.Dataset(options.FILE[0], 'r')

    var = f.variables['pism_config']

    print(header)

    for p in var.ncattrs():     # assume that this list is sorted
        if is_special(p):
            continue

        printers[var.getncattr(p + "_type")](var, p)
