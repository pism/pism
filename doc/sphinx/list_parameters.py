#!/usr/bin/env python

import netCDF4
import argparse


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
        return '%e' % number  # use scientific notation if a number is big
    elif int(number) == number:
        return '%d' % int(number)  # remove zeros after the decimal point
    elif abs(number) <= 1e-5:
        return '%e' % number  # use scientific notation if small (and not zero; previous case)
    else:
        return '%f' % number


def entry(name, T, value, option, description, choices=None):

    if len(value) == 0:
        value = "*no default*"

    if value.find(",") != -1:
        # A comma-separated list
        value = "``{}``".format(value.replace(",", ", "))

    if choices is not None:
        # This is always a comma-separated list.
        choices = "``{}``".format(choices.replace(",", ", "))
        # value is a keyword, so we should format it as a literal
        value = "``{}``".format(value)

    template = """
#. :config:`{name}` (*{T}*)

   :Value: {value}"""

    if choices is not None:
        template += """
   :Choices: {choices}"""

    if option is not None:
        option = ":opt:`{}`".format(option)
        template += """
   :Option: {option}"""

    template += """
   :Description: {description}"""

    return template.format(name=name,
                           T=T,
                           value=value,
                           choices=choices,
                           option=option,
                           description=description)


def print_string(var, name, T):
    print(entry(name,
                T,
                value(var, name),
                option(var, name),
                doc(var, name)))


def print_scalar(var, name, T):
    V = "{} ({})".format(number_to_string(value(var, name)), units(var, name))
    print(entry(name, T, V, option(var, name), doc(var, name)))


def print_integer(var, name, T):
    print(entry(name,
                T,
                str(value(var, name)),
                option(var, name),
                doc(var, name)))


def print_keyword(var, name, T):
    print(entry(name,
                T,
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


printers = {"string": print_string,
            "scalar": print_scalar,
            "integer": print_integer,
            "boolean": print_string,
            "keyword": print_keyword}

header = """.. -*- mode: rst -*-

.. DO NOT EDIT: This file was automatically generated. Edit src/pism_config.cdl instead.
"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("FILE", nargs=1)
    options = parser.parse_args()

    f = netCDF4.Dataset(options.FILE[0], 'r')

    var = f.variables['pism_config']

    print(header)

    for parameter in var.ncattrs():     # assume that this list is sorted
        if is_special(parameter):
            continue

        T = var.getncattr(parameter + "_type")

        printers[T](var, parameter, T)
