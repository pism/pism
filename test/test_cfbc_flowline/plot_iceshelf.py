#!/usr/bin/env python

"""This plot script for a flowline ice shelf setup visualizes
   the effect of the Calving Front Boundary Condition at ice shelves
   and of a modified driving stress scheme at the calving front.

   contact: reese@pik-potsdam.de and albrecht@pik-potsdam.de"""


#import MISMIP

from pylab import figure, subplot, hold, plot, xlabel, ylabel, title, axis, vlines, savefig, text, show, legend, xlim, ylim
from sys import exit

import numpy as np
from optparse import OptionParser
import os.path

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

# PISM defaults
rhoi=910.0
rhosw=1028
grav=9.81

ca=1750e3

def secpera():
    "Number of seconds per year."
    #return 3.15569259747e7
    #return 3600.0*24.0*365.0
    return 31536000.0


def process_options():
    "Process command-line options and arguments."
    parser = OptionParser()
    parser.usage = "%prog <input files> [options]"
    parser.description = "Plots profile along flowline."
    parser.add_option("-o", "--output", dest="output", type="string",
                      help="Output image file name (e.g. -o foo.png)")
    parser.add_option("-f", "--flux", dest="profile", action="store_false", default=True,
                      help="Plot ice flux")
    parser.add_option("-p", "--profile", dest="flux", action="store_false", default=True,
                      help="Plot geometry profile")
    parser.add_option("-d", "--taud", dest="taud", action="store_false", default=True,
                      help="Plot driving stress profile")


    opts, args = parser.parse_args()

    if len(args) == 0:
        print("ERROR: An input file is requied.")
        exit(0)

    if len(args) > 1 and opts.output:
        print("More than one input file given. Ignoring the -o option...\n")
        opts.output = None

    #if opts.output and opts.profile and opts.flux:
     #   print("Please choose between flux (-f) and profile (-p) plots.")
     #   exit(0)

    return args, opts.output, opts


def read(filename, name):
    "Read a variable and extract the middle row."
    nc = NC(filename)

    try:
        var = nc.variables[name][:]
    except:
        print("ERROR: Variable '%s' not present in '%s'" % (name, filename))
        exit(1)

    N = len(var.shape)
    if N == 1:
        return var[:]               # a coordinate variable ('x')
    elif N == 2:
        return var[1]               # get the middle row
    elif N == 3:
        return var[-1, 1]            # get the middle row of the last record
    else:
        raise Exception("Can't read %s. (It's %d-dimensional.)" % (name, N))


def plot_profile(in_file, out_file):
    #print("Reading %s to plot geometry profile")

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-profile.pdf"

    mask = read(in_file, 'mask')
    usurf = read(in_file, 'usurf')
    topg = read(in_file, 'topg')
    thk = read(in_file, 'thk')
    x = read(in_file, 'x')

    # mask out ice-free areas
    usurf = np.ma.array(usurf, mask=mask == 4)

    # compute the lower surface elevation
    lsrf = topg.copy()
    lsrf[mask == 3] = -rhoi / rhosw * thk[mask == 3]
    #lsrf[mask == 3] = -MISMIP.rho_i() / MISMIP.rho_w() * thk[mask == 3]
    lsrf = np.ma.array(lsrf, mask=mask == 4)

    # convert x to kilometers
    x /= 1e3

    figure(1)
    ax = subplot(111)
    hold(True)
    plot(x, np.zeros_like(x), ls='dotted', color='black')
    plot(x, topg, color='black')
    plot(x, usurf, '.-', color='grey', markersize=2,alpha=1.0)
    plot(x, lsrf,  '.-', color='grey', markersize=2,alpha=1.0)

    xlabel("x ($\mathrm{km}$)", size=14)
    ylabel(r"elevation ($\mathrm{m}$)", size=14)


    title("Flowline setup")

    _, _, ymin, ymax = axis(xmin=0, xmax=x.max())

    print("Saving '%s'...\n" % out_file)
    savefig(out_file)


def plot_flux(in_file, out_file):
    print("Reading %s to plot ice flux " % in_file) #for model %s, experiment %s, grid mode %s, step %s" % (
        #in_file, model, experiment, mode, step))

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-flux.pdf"

    x = read(in_file, 'x')
    flux_mag = read(in_file, 'flux_mag')

    # plot positive xs only
    flux_mag = flux_mag[x >= 0]
    x = x[x >= 0]

    figure(2)
    hold(True)

    plot(x / 1e3, flux_mag, 'k.-', markersize=10, linewidth=2)
    #plot(x / 1e3, x * MISMIP.a() * MISMIP.secpera(), 'r:', linewidth=1.5)

    #title("MISMIP experiment %s, step %d" % (experiment, step))
    xlabel("x ($\mathrm{km}$)", size=14)
    ylabel(r"flux ($\mathrm{m}^2\,\mathrm{a}^{-1}$)", size=14)

    print("Saving '%s'...\n" % out_file)
    savefig(out_file, dpi=300, facecolor='w', edgecolor='w')
    #show()


def plot_vel(dev_file, in_file, com_file, out_file):
    print("Reading %s to plot ice velocity" % (
        in_file)) #, model, experiment, mode, step))

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-ssavel.pdf"

    x = read(in_file, 'x')
    vel = read(in_file, 'u_ssa')
    thk = read(in_file, 'thk')

    xd = read(dev_file, 'x')
    veld = read(dev_file, 'u_ssa')
    thkd = read(dev_file, 'thk')

    x2 = read(comp_file, 'x')
    velc = read(comp_file, 'u_ssa_bc')

    # plot positive xs only
    veld = veld * secpera()
    vel = vel[x >= 0] * secpera()
    thk = thk[x >= 0]
    x = x[x >= 0] 

    figure(3)
    hold(True)

    plot(x / 1e3, veld, 'g.-', markersize=2, linewidth=1, alpha=0.5,label='dev')
    plot(xd / 1e3, vel, 'k.-', markersize=2, linewidth=1, alpha=0.5,label='fix')
    plot(x2 / 1e3, velc, 'r:', markersize=2, linewidth=1.5,label='veen')

    text(300,300,'max diff fix: '+str(np.max(vel)-np.max(velc))+' m/yr',color='k')
    text(300,600,'max diff dev: '+str(np.max(veld)-np.max(velc))+' m/yr',color='g')

    #title("MISMIP experiment %s, step %d" % (experiment, step))
    xlabel("x ($\mathrm{km}$)", size=14)
    ylabel(r"SSA vel ($\mathrm{m}\,\mathrm{a}^{-1}$)", size=14)

    legend()

    title("Velocity profile")

    print("Saving '%s'...\n" % out_file)
    savefig(out_file, dpi=300, facecolor='w', edgecolor='w')
    #show()

def plot_taud(dev_file, in_file, out_file):
    print("Reading %s to plot driving stress" % (
        in_file)) #, model, experiment, mode, step))

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-taud.pdf"

    x = read(in_file, 'x')
    taud = read(in_file, 'taud_mag')
    thk = read(in_file, 'thk')
    hs = read(in_file, 'usurf')
    td = np.zeros_like(taud)
    dx=x[1]-x[0]

    tauddev = read(dev_file, 'taud_mag')

    for i,xi in enumerate(x):
      if i>0 and i<len(x)-1:
        td[i]=-(hs[i+1]-hs[i-1])*thk[i]*rhoi*grav*0.5/dx


    figure(4)
    hold(True)

    plot(x / 1e3, taud, 'k.-', markersize=2, linewidth=1)
    plot(x / 1e3, tauddev, 'g.-', markersize=2, linewidth=1,alpha=0.5)
    plot(x / 1e3, td, 'r:', markersize=2, linewidth=1.5)

    xlim((ca-3*dx)/1e3,(ca+2*dx)/1e3)
    ylim(-5e3,5e3)

    #title("MISMIP experiment %s, step %d" % (experiment, step))
    xlabel("x ($\mathrm{km}$)", size=14)
    ylabel(r"SSA taud ($\mathrm{Pa}$)", size=14)

    print("Saving '%s'...\n" % out_file)
    savefig(out_file, dpi=300, facecolor='w', edgecolor='w')
    #show()





if __name__ == "__main__":
        args, out_file, opts = process_options()

        dev_file=args[0]
        in_file=args[1]
        comp_file=args[2]

        print dev_file,in_file,comp_file

        if opts.profile:
            plot_profile(in_file, out_file)

        #if opts.flux:
        #    plot_flux(in_file, out_file)

        if opts.taud:
            plot_taud(dev_file, in_file, out_file)

        
        plot_vel(dev_file, in_file, comp_file, out_file)

        #show()

