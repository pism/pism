src/
====

Do

        git ls-files | grep -P "(\.c|\.cc)" | xargs wc -l

to count lines of C and C++ code here.  As of Feb 2016, this gives about 72,000
lines, a fairly conservative count of PISM code lines because it excludes so
much test and example script code.

