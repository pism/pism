#!/usr/bin/env python
from sys import argv, exit
from getopt import getopt, GetoptError

def success():
    print "Files are the same within given tolerance."
    exit(0)

def failure():
    print "Files are different."
    exit(1)

def compare_vars(nc1, nc2, name, tol):
    from numpy import squeeze
    # Find variables:
    try:
        var1 = squeeze(nc1.variables[name][:])
        var2 = squeeze(nc2.variables[name][:])
    except:
        # This can happen if one of the files does not have the variable.
        failure()

    try:
        delta = abs(var1 - var2).max()
    except:
        # This can happen if variables have different shapes.
        failure()

    # The actual check:
    if (delta > tol):
        print "name = %s, delta = %e, tol = %e" % (name, delta, tol)
        failure()
    

def compare(file1, file2, variables, exclude, tol):
    try:
        from netCDF4 import Dataset as NC
    except:
        from netCDF3 import Dataset as NC

    from numpy import unique, r_

    try:
        nc1 = NC(file1, 'r')
        nc2 = NC(file2, 'r')
    except:
        # This can happen if one of the files could not be opened.
        failure()

    if (exclude == False):
        if len(variables) == 0:
            vars1 = nc1.variables.keys()
            vars2 = nc2.variables.keys()
            variables = unique(r_[vars1, vars2])

        for each in variables:
            compare_vars(nc1, nc2, each, tol)
    else:
        vars1 = nc1.variables.keys()
        vars2 = nc2.variables.keys()
        vars = unique(r_[vars1, vars2])

        for each in vars:
            if (each in variables):
                continue
            compare_vars(nc1, nc2, each, tol)


if __name__ == "__main__":
    from numpy import double
    opts, args = getopt(argv[1:], "t:v:x")
    tol = 0
    file1 = ""
    file2 = ""
    variables = []
    exclude = False
    for (opt, arg) in opts:
        if opt == "-t":
            tol = double(arg)
        if opt == "-x":
            exclude = True
        if opt == "-v":
            variables = arg.split(",")

    if len(args) != 2:
        failure()

    compare(args[0],args[1], variables, exclude, tol)

    success()

    
            
