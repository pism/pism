#! /usr/bin/env python
## This script takes the standard out from verifynow.py and produces graphs using gnuplot.

import getopt
import sys
import numpy
import Gnuplot

# defaults
IN_FILE = 'foo.txt'
OUT_FILE_START = 'report'
TEST_NAME = 'A'
EXCLUDE = ''

##### process command line arguments #####
infilename = IN_FILE
outfilestart = OUT_FILE_START
testname = TEST_NAME
errtype = 'geometry'
exclude = EXCLUDE
try:
  opts, args = getopt.getopt(sys.argv[1:], "f:t:o:e:x:",
                             ["file=", "test=", "out=", "error=", "exclude="])
  for opt, arg in opts:
    if opt in ("-f", "--file"):
      infilename = arg
    elif opt in ("-t", "--test"):
      testname = arg
    elif opt in ("-o", "--out"):
      outfilestart = arg
    elif opt in ("-e", "--error"):
      errtype = arg
    elif opt in ("-x", "--exclude"):
      exclude = arg
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
  print 'FAILED: no errors of type "' + errtype + '" found'
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

# plot to postscript (.eps) files
g = Gnuplot.Gnuplot()
g('set data style linespoints')
g('set pointsize 2')
g('set logscale xy 10')
g('set key top left')
g.xlabel(dxname + '  (' + units + ')')
for nn in range(len(tags)):
  if exclude.find(tags[nn]) < 0: # if the tag is not in the exclude string, then plot
    d = Gnuplot.Data( path[0:len(err[nn])], err[nn], title=(tags[nn] + ' error') )
    if nn == 0:
      g.plot(d)
    else:
      g.replot(d)
  else:
    print 'excluding "' + tags[nn] + '" errors from plot'
g.title(titles[0] + ' numerical errors in test ' + testname)
outfilename = outfilestart + '_' + testname + '.eps'
g.hardcopy(outfilename, mode='eps', enhanced=1, color=1)
print 'postscript file ' + outfilename + ' written'
# done

