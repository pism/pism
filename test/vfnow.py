#!/usr/bin/env python

# @package vfnow
# \author Ed Bueler and Constantine Khroulev, University of Alaska Fairbanks, USA
# \brief A script for verification of numerical schemes in PISM.
# \details It specifies a refinement path for each of Tests ABCDEFGIJKL and runs
# pismv accordingly.
# Copyright (C) 2007--2013, 2015, 2016, 2018 Ed Bueler and Constantine Khroulev
##
# Organizes the process of verifying PISM.  It specifies standard refinement paths for each of the tests described in the user manual.  It runs the tests, times them, and summarizes the numerical errors reported at the end.
##
# Examples:
# - \verbatim vfnow.py \endverbatim use one processor and do three levels of refinement; this command is equivalent to \verbatim vfnow.py -n 2 -l 2 -t CGIJ \endverbatim,
# - \verbatim vfnow.py -n 8 -l 5 -t J --prefix=bin/ --mpido='aprun -n' \endverbatim will use \verbatim aprun -n 8 bin/pismv \endverbatim as the command and do five levels (the maximum) of refinement only on test J,
# - \verbatim vfnow.py -n 2 -l 3 -t CEIJGKLO \endverbatim uses two processers (cores) and runs in about an hour,
# - \verbatim vfnow.py -n 40 -l 5 -t ABCDEFGIJKLO \endverbatim will use forty processors to do all possible verification as managed by \c vfnow.py; don't run this unless you have a big computer and you are prepared to wait.
# For a list of options do \verbatim test/vfnow.py --help \endverbatim.
# Timing information is given in the \c vfnow.py output so performance, including parallel performance, can be assessed along with accuracy.

import sys
import time
import subprocess
from numpy import array

# A class describing a refinement path and command-line options
# for a particular PISM verification test.


class PISMVerificationTest:

    # max number of levels that will work with
    N = 50
    # one-letter test name
    name = ""
    # test description
    test = ""
    # description of the refinement path
    path = ""
    Mx = []
    My = []
    # 31 levels in the ice
    Mz = [31] * N
    # no bedrock by default
    Mbz = [1] * N
    # extra options (such as -y, -ys, -ssa_rtol)
    opts = ""
    executable = "pismv"

    def build_command(self, exec_prefix, level):
        M = list(zip(self.Mx, self.My, self.Mz, self.Mbz))

        if level > len(M):
            print("Test %s: Invalid refinement level: %d (only %d are available)" % (
                self.name, level, len(M)))
            return ""

        grid_options = "-Mx %d -My %d -Mz %d -Mbz %d" % M[level - 1]
        return "%s%s -test %s %s %s" % (exec_prefix, self.executable, self.name, grid_options, self.opts)


def run_test(executable, name, level, extra_options="", debug=False):
    try:
        test = tests[name]
    except:
        print("Test %s not found." % name)
        return

    if level == 1:
        print("# ++++  TEST %s:  verifying with %s exact solution  ++++\n# %s" % (
            test.name, test.test, test.path))
    else:
        extra_options += " -append"

    command = test.build_command(executable, level) + " " + extra_options

    if debug:
        print('# L%d\n%s' % (level, command))
        return
    else:
        print(' L%d: trying "%s"' % (level, command))

    # run PISM:
    try:
        lasttime = time.time()
        (status, output) = subprocess.getstatusoutput(command)
        elapsetime = time.time() - lasttime
    except KeyboardInterrupt:
        sys.exit(2)
    if status:
        sys.exit(status)
    print(' finished in %7.4f seconds; reported numerical errors as follows:' % elapsetime)

    # process the output:
    position = output.find('NUMERICAL ERRORS')
    if position >= 0:
        report = output[position:output.find('NUM ERRORS DONE')]
        endline = report.find('\n')
        print('    ' + report[0:endline])
        report = report[endline + 1:]
        while (len(report) > 1) and (endline > 0):
            endline = report.find('\n')
            if endline == -1:
                endline = len(report)
            print('   #' + report[0:endline])
            report = report[endline + 1:]
            endline = report.find('\n')
            if endline == -1:
                endline = len(report)
            print('   |' + report[0:endline])
            report = report[endline + 1:]
    else:
        print(" ERROR: can't find reported numerical error")
        sys.exit(99)


def define_refinement_paths(KSPRTOL, SSARTOL):
    # Define all the supported refinement paths:
    tests = {}
    # A
    A = PISMVerificationTest()
    A.name = "A"
    A.test = "steady, marine margin isothermal SIA"
    A.path = "(refine dx=53.33,40,26.67,20,13.33,km, dx=dy and Mx=My=31,41,61,81,121)"
    A.Mx = [31, 41, 61, 81, 121]
    A.My = A.Mx
    A.opts = "-y 25000.0"
    tests['A'] = A
    # B
    B = PISMVerificationTest()
    B.name = "B"
    B.test = "moving margin isothermal SIA (Halfar)"
    B.path = "(refine dx=80,60,40,30,20,km, dx=dy and Mx=My=31,41,61,81,121)"
    B.Mx = [31, 41, 61, 81, 121]
    B.My = B.Mx
    B.opts = "-ys 422.45 -y 25000.0"
    tests['B'] = B
    # C
    C = PISMVerificationTest()
    C.name = "C"
    C.test = "non-zero accumulation moving margin isothermal SIA"
    C.path = "(refine dx=50,33.33,25,20,16,km, dx=dy and Mx=My=41,61,81,101,121)"
    C.Mx = [41, 61, 81, 101, 121]
    C.My = C.Mx
    C.opts = "-y 15208.0"
    tests['C'] = C
    # D
    D = PISMVerificationTest()
    D.name = "D"
    D.test = "time-dependent isothermal SIA"
    D.path = "(refine dx=50,33.33,25,20,16.67,km, dx=dy and Mx=My=41,61,81,101,121)"
    D.Mx = [41, 61, 81, 101, 121]
    D.My = D.Mx
    D.opts = "-y 25000.0"
    tests['D'] = D
    # E
    E = PISMVerificationTest()
    E.name = "E"
    E.test = "steady sliding marine margin isothermal SIA"
    E.path = "(refine dx=53.33,40,26.67,20,13.33,km, dx=dy and Mx=My=31,41,61,81,121)"
    E.Mx = [31, 41, 61, 81, 121]
    E.My = E.Mx
    E.opts = "-y 25000.0"
    tests['E'] = E
    # F
    F = PISMVerificationTest()
    F.name = "F"
    F.test = "steady thermomechanically-coupled SIA"
    F.path = "(refine dx=30,20,15,10,7.5,km, dx=dy, dz=66.67,44.44,33.33,22.22,16.67 m and Mx=My=Mz=61,91,121,181,241)"
    F.Mx = [61, 91, 121, 181, 241]
    F.My = F.Mx
    F.Mz = F.Mx
    F.opts = "-y 25000.0"
    tests['F'] = F
    # G
    G = PISMVerificationTest()
    G.name = "G"
    G.test = "time-dependent thermomechanically-coupled SIA"
    G.path = "(refine dx=30,20,15,10,7.5,km, dx=dy, dz=66.67,44.44,33.33,22.22,16.67 m and Mx=My=Mz=61,91,121,181,241)"
    G.Mx = [61, 91, 121, 181, 241]
    G.My = G.Mx
    G.Mz = G.Mx
    G.opts = "-y 25000.0"
    tests['G'] = G
    # H
    H = PISMVerificationTest()
    H.name = "H"
    H.test = "moving margin, isostatic bed, isothermal SIA"
    H.path = "(refine dx=80,60,40,30,20,km, dx=dy and Mx=My=31,41,61,81,121)"
    H.Mx = [31, 41, 61, 81, 121]
    H.My = H.Mx
    H.opts = "-bed_def iso -y 60000.0"
    tests['H'] = H
    # I
    I = PISMVerificationTest()
    I.executable = "ssa_testi"
    I.name = "I"
    I.test = "plastic till ice stream (SSA)"
    I.path = "(refine dy=5000,1250,312.5,78.13,19.53,m, My=49,193,769,3073,12289)"
    I.Mx = [5] * 5
    I.My = [49, 193, 769, 3073, 12289]
    I.executable = "ssa_testi"
    I.opts = "-ssa_method fd -ssa_rtol %1.e -ssafd_ksp_rtol %1.e" % (SSARTOL, KSPRTOL)
    tests['I'] = I
    # J
    J = PISMVerificationTest()
    J.executable = "ssa_testj"
    J.name = "J"
    J.test = "periodic ice shelf (linearized SSA)"
    J.Mx = [49, 98, 196, 392, 784]
    J.My = J.Mx
    J.path = "(refine Mx={})".format(J.Mx)
    J.Mz = [11] * 5
    J.executable = "ssa_testj"
    J.opts = "-ssa_method fd -ssafd_pc_type asm -ssafd_sub_pc_type lu -ssafd_ksp_rtol %1.e" % KSPRTOL
    tests['J'] = J
    # K
    K = PISMVerificationTest()
    K.name = "K"
    K.test = "pure conduction problem in ice and bedrock"
    K.path = "(refine dz=100,50,25,12.5,6.25,m, Mz=41,81,161,321,641)"
    K.Mx = [8] * 5
    K.My = K.Mx
    K.Mz = array([41, 81, 161, 321, 641])
    K.Mbz = (K.Mz - 1) / 4 + 1
    K.opts = "-y 130000.0 -Lbz 1000 -z_spacing equal"
    tests['K'] = K
    # L
    L = PISMVerificationTest()
    L.name = "L"
    L.test = "non-flat bed stead isothermal SIA"
    L.path = "(refine dx=60,30,20,15,10,km, dx=dy and Mx=My=31,61,91,121,181)"
    L.Mx = [31, 61, 91, 121, 181]
    L.My = L.Mx
    L.opts = "-y 25000.0"
    tests['L'] = L
    # M
    M = PISMVerificationTest()
    M.name = "M"
    M.test = "annular ice shelf with a calving front (SSA)"
    M.path = "(refine dx=50,25,16.666,12.5,8.333 km; dx=dy and My=31,61,91,121,181)"
    M.Mx = [31, 61, 91, 121, 181]
    M.My = M.Mx
    M.Mz = [11] * 5
    M.opts = "-ssa_rtol %1.e -ssafd_ksp_rtol %1.e" % (SSARTOL, KSPRTOL)
    tests['M'] = M
    # O
    O = PISMVerificationTest()
    O.name = "O"
    O.test = "basal melt rate from conduction problem in ice and bedrock"
    O.path = "(refine dz=100,50,25,12.5,6.25,m, Mz=41,81,161,321,641)"
    O.Mx = [8] * 5
    O.My = O.Mx
    O.Mz = array([41, 81, 161, 321, 641])
    O.Mbz = (O.Mz - 1) / 4 + 1
    O.opts = "-z_spacing equal -zb_spacing equal -Lbz 1000 -y 1000 -no_mass"
    tests['O'] = O

    # test K (for a figure in the User's Manual)
    K = PISMVerificationTest()
    K.name = "K"
    K.test = "pure conduction problem in ice and bedrock"
    K.path = "(lots of levels)"
    K.Mz = array([101, 121, 141, 161, 181, 201, 221, 241, 261, 281, 301, 321])
    K.Mbz = (K.Mz - 1) / 4 + 1
    K.Mx = [8] * len(K.Mz)
    K.My = K.Mx
    K.opts = "-y 130000.0 -Lbz 1000"
    tests['K_manual'] = K

    # test B (for a figure in the User's Manual)
    B = PISMVerificationTest()
    B.name = "B"
    B.test = "moving margin isothermal SIA (Halfar)"
    B.path = "(lots of levels)"
    B.Mx = [31, 41, 51, 61, 71, 81, 91, 101, 111, 121]
    B.My = B.Mx
    B.Mz = [31] * len(B.Mx)
    B.Mbz = [1] * len(B.Mx)
    B.opts = "-ys 422.45 -y 25000.0"
    tests['B_manual'] = B

    # test G (for a figure in the User's Manual)
    G = PISMVerificationTest()
    G.name = "G"
    G.test = "time-dependent thermomechanically-coupled SIA"
    G.path = "(lots of levels)"
    G.Mx = [61, 71, 81, 91, 101, 111, 121, 151, 181]
    G.My = G.Mx
    G.Mz = G.Mx
    G.opts = "-y 25000.0"
    tests['G_manual'] = G

    # test I (for a figure in the User's Manual)
    I = PISMVerificationTest()
    I.executable = "ssa_testi"
    I.name = "I"
    I.test = "plastic till ice stream (SSA)"
    I.path = "(lots of levels)"
    I.My = [51, 101, 151, 201, 401, 601, 801, 1001, 1501, 2001, 2501, 3073]
    I.Mx = [5] * len(I.My)
    I.opts = "-ssa_method fd -ssa_rtol %1.e -ssafd_ksp_rtol %1.e" % (SSARTOL, KSPRTOL)
    tests['I_manual'] = I

    return tests


from argparse import ArgumentParser

parser = ArgumentParser()
parser.description = """PISM verification script"""
parser.add_argument("--eta", dest="eta", action="store_true",
                    help="to add '-eta' option to pismv call")
parser.add_argument("-l", dest="levels", type=int, default=2,
                    help="number of levels of verification; '-l 1' fast, '-l 5' slowest")
parser.add_argument("--mpido", dest="mpido", default="mpiexec -np",
                    help="specify MPI executable (e.g. 'mpirun -np' or 'aprun -n')")
parser.add_argument("-n", dest="n", type=int, default=2,
                    help="number of processors to use")
parser.add_argument("--prefix", dest="prefix", default="",
                    help="path prefix to pismv executable")
parser.add_argument("-r", dest="report_file", default="",
                    help="name of the NetCDF error report file")
parser.add_argument("-t", dest="tests", nargs="+",
                    help="verification tests to use (A,B,C,D,E,F,G,H,I,J,K,L,M,O); specify a space-separated list", default=['C', 'G', 'I', 'J'])
parser.add_argument("-u", dest="unequal", action="store_true",
                    help="use quadratic vertical grid spacing")
parser.add_argument("--debug", dest="debug", action="store_true",
                    help="just print commands in sequence (do not run pismv)")
parser.add_argument("--manual", dest="manual", action="store_true",
                    help="run tests necessary to produce figures in the User's Manual")

options = parser.parse_args()
extra_options = ""

if options.eta:
    extra_options += " -eta"

if options.unequal:
    extra_options += " -z_spacing quadratic"

if options.report_file:
    extra_options += " -report_file %s" % options.report_file

predo = ""
if options.n > 1:
    predo = "%s %d " % (options.mpido, options.n)

exec_prefix = predo + options.prefix

KSPRTOL = 1e-12  # for tests I, J, M
SSARTOL = 5e-7  # ditto
tests = define_refinement_paths(KSPRTOL, SSARTOL)

manual_tests = ["B_manual", "G_manual", "K_manual", "I_manual"]
if options.manual:
    print("# VFNOW.PY: test(s) %s, using '%s...'\n" % (manual_tests, exec_prefix) +
          "#           and ignoring options -t and -l")
    for test in manual_tests:
        N = len(tests[test].Mx)
        for j in range(1, N + 1):
            run_test(exec_prefix, test, j, extra_options,
                     options.debug)
else:
    print("# VFNOW.PY: test(s) %s, %d refinement level(s), using '%s...'" % (
        options.tests, options.levels, exec_prefix))

    for test in options.tests:
        for j in range(1, options.levels + 1):
            run_test(exec_prefix, test, j, extra_options,
                     options.debug)
