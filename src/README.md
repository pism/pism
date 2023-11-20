src/
====

Do

        git ls-files | grep -P "(\.c|\.cc)" | xargs wc -l

to count lines of C and C++ code here.  As of October 2017, this gives about
70,000 lines, a conservative count of PISM code lines because it excludes all
test/example script codes.
