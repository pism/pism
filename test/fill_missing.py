#! /usr/bin/env python

from pycdf import *
from numpy import *
import sys
import commands
import getopt

def get_bad(a, bad_value):
  """Scans the matrix 'a' and returns the indicies where 'bad_value'
  occurs."""
  points = []
  s = shape(a)

  for i in range(s[0]):
    for j in range(s[1]):
      if a[i][j] == bad_value:
        points.append(array((i, j)))

  return points

def laplace(a, bad_value, eps=1.):
  """Computes the solution to laplace's equation in the region where the
  values of a are equal to bad_value. Boundary conditions around bad_value
  regions are Dirichlet."""
  t = transpose
  ind = get_bad(a, bad_value)

  change = 100
  #eps = 1

  dimensions = shape(a)

  #before = sqrt((a[tuple(t(ind))]**2).mean())
  before = sqrt( (a[tuple(t(ind))]**2).mean() )
  while (abs(change) > eps):
    for e in ind:
      e = tuple(e)
      row = [1, -1 , 0, 0] + e[0]
      col = [0, 0, 1, -1] + e[1]
      
      # this section gets rid of the indices that are
      # out of bounds
      cond = greater_equal(row, 0)
      row = tuple(t(compress(cond, row)))
      col = tuple(t(compress(cond, col)))

      cond = less(row, dimensions[0])
      row = tuple(t(compress(cond, row)))
      col = tuple(t(compress(cond, col)))

      cond = greater_equal(col, 0)
      col = tuple(t(compress(cond, col)))
      row = tuple(t(compress(cond, row)))

      cond = less(col, dimensions[1])
      col = tuple(t(compress(cond, col)))
      row = tuple(t(compress(cond, row)))
     
      row = tuple(row)
      col = tuple(col)
      u = (row, col)
      a[tuple(e)] = a[u].mean()
    after = sqrt( (a[tuple(t(ind))]**2).mean() )
    #after = sqrt((a[tuple(t(ind))]**2).mean())
    change = after-before
    before = after
    print "change: " + str(change)
        
  return size(ind)

# command line arguments set
fileSet = 0
variablesSet = 0
outFileSet = 0
fileName = "guess.nc"
outFileName = "fill_missing_out.nc"
variables = ("bed")

try:
  opts, args = getopt.getopt(sys.argv[1:], "i:v:o:", ["in_file", "variables", "out_file"])
  for opt, arg in opts:
    if opt in ("-i", "--in_file"):
      fileName = arg
      fileSet = 1;
    if opt in ("-o", "--out_file"):
      outFileName = arg
      outFileSet = 1
    if opt in ("-v", "--variables"):
      variables = arg.split(",");
      variablesSet = 1
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

if fileSet == 0:
  print 'ERROR: No file name given, exiting...'
  sys.exit(2)
if variablesSet == 0:
  print 'ERROR: No variable names given, exiting...'
  sys.exit(2)

try:
  (status, output)=commands.getstatusoutput('cp '+fileName+' '+outFileName)
except:
  print 'output: '+output
  sys.exit(status)

print 'Reading from file: ' + fileName

try:
  ncfile = CDF(outFileName, NC.WRITE)
  ncfile.automode()
except CDFError,err:
  print err[2]
  sys.exit(2)

for varName in variables:
  try:
    var = ncfile.var(varName)
    vMissing = var.attr("missing_value")
    bad_value = vMissing.get(NC.FLOAT)
    print "Smoothing "+ varName + " with missing value: " + str(bad_value)
    b = var[:]
    num = laplace(b, float(bad_value))
    var[:] = b
    print "Number of missing values for "+ varName + ": " + str(num)
  except CDFError,err:
    print "ERROR, variable " + varName + " reported the following error: \""\
          + err[2] + "\""

ncfile.close()

