#!/usr/bin/env python
import os
import sys

print("Test # 0: presence of tools and Python modules needed by other tests.")

tools = ["ncpdq", "ncap2", "ncks", "ncrename", "diff", "rm", "cp"]
for tool in tools:
    print("Looking for %s..." % tool)
    if os.system("which %s" % tool) != 0:
        print("ERROR: %s does not seem to be on the PATH!" % tool)
        sys.exit(1)

modules = ["numpy", "sys", "netCDF4"]
for module in modules:
    try:
        print("Trying to import %s..." % module)
        __import__(module)
        print("OK.")
    except:
        print("\nERROR: Python cannot import %s" % module)
        sys.exit(1)
