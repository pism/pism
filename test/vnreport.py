#! /usr/bin/env python
## This script takes the standard out from verifynow.py and produces graphs.

import getopt
import sys
import numpy
import Gnuplot

# defaults
IN_FILE = 'foo.txt'
OUT_FILE_START = 'report'
TEST_NAME = 'A'

##### process command line arguments #####
pref = ''
infilename = IN_FILE
outfilestart = OUT_FILE_START
testname = TEST_NAME
errtype = 'geometry'
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:f:t:o:e:",
                             ["prefix=", "file=", "test=", "out=", "error="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      pref = arg
    elif opt in ("-f", "--file"):
      infilename = arg
    elif opt in ("-t", "--test"):
      testname = arg
    elif opt in ("-o", "--out"):
      outfilestart = arg
    elif opt in ("-e", "--error"):
      errtype = arg
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

# read file
print "reading from " + infilename + " (which was standard output from verifynow.py, right?)"
try:
  infile = open(infilename, 'r');
except IOError:
  print 'ERROR: File: ' + infilename + ' could not be found.'
  sys.exit(2)

# search for numerical errors from chosen test
while True:
  myline = infile.readline()
  if not myline: # stop if nothing left to read
    foundmytest = False
    break
  if myline.find('**TEST ' + testname) >= 0:
    foundmytest = True
    break
if foundmytest == False:
  print 'FAILED: numerical errors for test ' + testname + ' not found in file ' + infilename
  sys.exit(3)

# the next line will contain the refinement path
myline = infile.readline()
pos = myline.find('path ')
if pos < 0:
  print 'FAILED: refinement path not found in expected spot'
  sys.exit(5)
dxname = myline[pos+5:pos+7]
myline = myline[pos+8:]
strpath = myline.split(',')  # split on COMMAS!
path = []
for kk in range(5):  # expect exactly 5 pts on refinement path!
  path.append(float(strpath[kk]))
units = strpath[5]  # expect units to be given next

# read and parse error report into tokens
tags=[]   # the names of the error values
vals=[]   # the values themselves
titles=[] # category of the error; will be title on plot
count = 0
while True:
  myline = infile.readline()
  # stop if nothing left to read or we run into another test report
  if (not myline) or (myline.find('**TEST') >= 0): 
    break
  while myline[0] == ' ': # strip off leading spaces
    myline = myline[1:]
  if myline[0] == '#':
    titles.append(myline[1:11])
    myline = myline[12:] #strip off 12
    tags.append(myline.split())
  elif myline[0] == '|':
    myline = myline[12:] #strip off 12
    myvals = []
    for num in myline.split():
      myvals.append(float(num))
    vals.append(myvals)
    count = count + len(myvals)
print str(count) + ' total numerical error values read for test ' + testname
infile.close()

# go through and eliminate all but the desired error type
NN = 0
while len(titles) > NN:
  if titles[NN].find(errtype) >= 0:
    NN = NN + 1
  else:
    #remove title and corresponding tags and errors
    titles = titles[0:NN] + titles[NN+1:]
    tags = tags[0:NN] + tags[NN+1:]
    vals = vals[0:NN] + vals[NN+1:]
if NN == 0:
  print 'FAILED: no errors of given type found'
  sys.exit(5)

# check that tags are consistent; then discard redundant tags
for kk in range(len(tags[0])):
  firsttag = tags[0][kk]
  for ll in range(len(tags)):
    if not (firsttag == tags[ll][kk]):
      print 'FAILED: inconsistent tags'
      sys.exit(4)
tags = tags[0]

# reorganize errors; now err[N] is along refinement path
err = numpy.array(vals).transpose()
print 'found ' + str(len(err)*len(err[0])) + \
      ' "' + errtype + '"-type  numerical error values for plotting'

# plot to postscript files
g = Gnuplot.Gnuplot()
g.title('numerical errors in test ' + testname + ' ' + titles[0])
g('set data style linespoints')
g('set logscale xy 10')
g.xlabel(dxname + '  (' + units + ')')
for nn in range(len(tags)):
  d = Gnuplot.Data(path[0:len(err[nn])], err[nn])
  g.ylabel(tags[nn] + ' error')
  g.plot(d)
  outfilename = outfilestart + '_' + testname + '_' + tags[nn] + '.ps'
  g.hardcopy(outfilename, enhanced=1, color=1)
  print 'postscript file ' + outfilename + ' written'
# done

