#!/usr/bin/env python
from sys import argv, exit, stdin
from getopt import getopt, GetoptError

usage="""
BELOW.PY  Print all lines from a text file in which entries in a particular
column are numbers below a given level.

Stops at the end of the file OR the first line which contains only whitespace.

Examples:

1) To show lines in diffs.txt where column 4 contains values less than
0.5, do:

  ./below.py -c 4 -l 0.5 diffs.txt

2) To show lines in diffs.txt where column 2 contains values which are
less than 0.01 in absolute value, do:

  ./below.py -a -c 2 -l 0.01 diffs.txt
  
3) To count the lines:

  ./below.py -a -c 2 -l 0.01 diffs.txt | wc -l
  
4) Chain together.  To show lines in diffs.txt where column 2 is less that 0.01
in absolute value AND column 3 is less than 1.0:

  ./below.py -a -c 2 -l 0.01 diffs.txt | ./below.py -c 4 -l 1.0
"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

if __name__ == "__main__":
    try:
      opts, args = getopt(argv[1:], "ac:l:", ["help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR below.py')
    useabsolute = False
    level = 0.0
    col = 1
    for (opt, optarg) in opts:
        if opt == "-a":
            useabsolute = True
        if opt == "-c":
            col = int(optarg)
        if opt == "-l":
            level = float(optarg)
        if opt in ("--help", "--usage"):
            print usage
            exit(0)
    #print "col = %d, level = %f" % (col,level)
    
    if len(args) == 0:
      ff = stdin
      #print "  reading lines from stdin"
    elif len(args) == 1:
      try:
        ff = file(args[0],'r')
      except:
        usagefailure('ERROR: could not read from file %s' % args[0])
      #print "  reading lines from text file %s and ..." % args[0]
    else:
      usagefailure('ERROR: below.py requires exactly one file name, or it reads from stdin')

    if col <= 0:
      usagefailure('ERROR: option -c must give column number at least one')

    ffline = ff.readline()
    fflen = len(ffline.split())
    while fflen > 0:
      if fflen < col:
        print "ERROR:  found a line with less than %d columns" % col
        exit(1)
      entry = ffline.split()[col-1]
      if useabsolute:
        if abs(float(entry)) < level:
          print ffline,
      else:
        if float(entry) < level:
          print ffline,
      ffline = ff.readline()
      fflen = len(ffline.split())

    if len(args) != 0:
      ff.close()

