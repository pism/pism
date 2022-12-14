#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as NC
from scipy.integrate import quadrature

import bedrock_step as exact

def V_exact():
    "Exact volume"
    def H(x):
        return exact.surface(x) - exact.bed_elevation(x)

    # Split the interval [0, x_m] into [0, x_s] and [x_s, x_m] because H(x) has a
    # discontinuity at x_s:
    V1, err1 = quadrature(H, 0, exact.x_s, maxiter=100)
    V2, err2 = quadrature(H, exact.x_s, exact.x_m, maxiter=100)

    return V1 + V2, err1 + err1

def plot(filename, output):

    with NC.Dataset(filename, "r") as f:
        x = f.variables["x"][:]
        # dimensions: time, y, x
        s_model = f.variables["usurf"][-1, 1, :]

        subset = np.logical_and(x >= 0, x <= 3e4)

        s_model = s_model[subset]
        x = x[subset]

    bed_exact = exact.bed_elevation(x)
    s_exact = exact.surface(x)

    fig, ax = plt.subplots()

    V_model = np.trapz(s_model - bed_exact, x)
    V_e, err = V_exact()

    print(f"Exact volume: {V_e} m^2 (error: {err} m^2)")
    print(f"Modeled volume: {V_model} m^2")
    print(f"Relative volume error: {(V_model - V_e) / V_e * 100} %")

    ax.plot(x/1000, bed_exact, label="bed", color="black")
    ax.plot(x/1000, s_exact, "--", color="gray", label="exact steady state surface elevation")
    ax.plot(x/1000, s_model, label="modeled surface elevation")
    ax.set_xlabel('x [km]')
    ax.set_ylabel('z [m]')
    ax.set_title("Modeled ice geometry at 50000 years")
    ax.legend()

    fig.set_size_inches(8, 4)
    fig.savefig(output)

if __name__ == "__main__":
    import sys
    plot(sys.argv[1], sys.argv[2])
