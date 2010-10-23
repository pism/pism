#!/usr/bin/env python
from sys import argv, exit, stdin
from getopt import getopt, GetoptError
from numpy import array

usage="""
WEIGHTS.PY  Read a text file with real numbers in columns 2,3,4.  Computes a 
weighted sum of those numbers with the given weights.  Gives special 
interpretation to column 2: first the value is squared (and then weighted).

Example:
  ./weights.py -2 1.0 -3 1.0 -4 1.0 diffs.txt Wdiffs.txt

Run with --help or --usage to get a this usage message.
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

if __name__ == "__main__":
    try:
      opts, args = getopt(argv[1:], "2:3:4:", ["help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR objective.py')
    w = array([1.0, 1.0, 1.0])
    for (opt, optarg) in opts:
        if opt == "-2":
            w[0] = float(optarg)
        if opt == "-3":
            w[1] = float(optarg)
        if opt == "-4":
            w[2] = float(optarg)
        if opt in ("--help", "--usage"):
            print usage
            exit(0)
    
    if len(args) != 2:
      usagefailure('ERROR: weights.py requires exactly two file names')

    ff = file(args[0], 'r')
    fout = file(args[1], 'a')
    
    while True:
      ffline = ff.readline()
      fflen = len(ffline.split())
      if fflen <= 0:
        break
      if fflen < 4:
        print "ERROR:  found a line in input file with fewer than 4 columns"
        print "EXITING EARLY ..."

      name, col2, col3, col4 = ffline.split()
      wsum = w[0] * float(col2)**2.0 + w[1] * float(col3) + w[2] * float(col4)
      
      fout.write("%s %12.7f\n" % (ffline[:-1], wsum))

    ff.close()
    fout.close()

