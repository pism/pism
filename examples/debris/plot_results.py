#!/usr/bin/env python3
"""Plot the debris-covered valley glacier example: maps of the glacier and the debris
cover, longitudinal profiles along the valley center line, and the debris mass budget.

Usage: plot_results.py spatial.nc scalar.nc figure.png
"""

import sys
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def last(var):
    return np.squeeze(var[-1, :, :])


def main(spatial_file, scalar_file, figure):
    with netCDF4.Dataset(spatial_file) as f:
        x = f.variables["x"][:] / 1e3
        y = f.variables["y"][:] / 1e3
        thk = last(f.variables["thk"])
        usurf = last(f.variables["usurf"])
        speed = last(f.variables["velsurf_mag"])
        h_d = last(f.variables["debris_thickness"])
        cover = last(f.variables["debris_cover_fraction"])
        column = last(f.variables["englacial_debris_column_mass"])
        C = f.variables["englacial_debris_concentration"][-1, :, :, :]
        z = f.variables["z"][:]
        time = f.variables["time"][-1]

    with netCDF4.Dataset(scalar_file) as f:
        t = np.asarray(f.variables["time"][:], dtype=float)
        units = f.variables["time"].units
        if units.startswith("seconds"):
            t = t / (365.0 * 86400.0)
        elif units.startswith("days"):
            t = t / 365.0
        M_eng = f.variables["englacial_debris_mass"][:]
        M_sup = f.variables["supraglacial_debris_mass"][:]
        M_tot = f.variables["total_debris_mass"][:]
        # the fluxes are written in "glaciological" units (kg / year)
        def per_year(name):
            v = f.variables[name]
            return v[:] * (1.0 if "year" in v.units else 365.0 * 86400.0)
        F_in = per_year("debris_input_mass_flux")
        F_melt = per_year("debris_melt_out_mass_flux")
        F_out = per_year("debris_output_mass_flux")
        err = f.variables["debris_mass_conservation_error"][:]

    j_c = len(y) // 2
    icy = thk > 1.0

    fig, axes = plt.subplots(3, 2, figsize=(12, 11))

    ax = axes[0, 0]
    m = ax.pcolormesh(x, y, np.ma.array(speed, mask=~icy), shading="auto")
    ax.contour(x, y, thk, levels=[1.0], colors="k", linewidths=0.5)
    fig.colorbar(m, ax=ax, label="surface speed (m/yr)")
    ax.set_title(f"glacier after {float(time) / (365 * 86400):.0f} years")
    ax.set_aspect("equal")

    ax = axes[0, 1]
    m = ax.pcolormesh(x, y, np.ma.array(h_d, mask=~icy), shading="auto", cmap="YlOrBr")
    ax.contour(x, y, thk, levels=[1.0], colors="k", linewidths=0.5)
    fig.colorbar(m, ax=ax, label="debris thickness (m)")
    ax.set_title("supraglacial debris")
    ax.set_aspect("equal")

    ax = axes[1, 0]
    ax.plot(x, usurf[j_c, :], label="surface")
    ax.plot(x, usurf[j_c, :] - thk[j_c, :], "k", label="bed")
    ax.set_ylabel("elevation (m)")
    ax.set_xlabel("x (km)")
    ax.legend()
    ax2 = ax.twinx()
    ax2.plot(x, h_d[j_c, :], "C3")
    ax2.set_ylabel("debris thickness (m)", color="C3")

    ax = axes[1, 1]
    # englacial concentration along the center line, as a function of height above the bed
    # (PISM stores 3D fields as (y, x, z))
    profile = np.asarray(C[j_c, :, :]).T        # (z, x)
    Hc = thk[j_c, :]
    Z = z[:, None] * np.ones_like(Hc)[None, :]
    profile = np.ma.array(profile, mask=Z > Hc[None, :])
    m = ax.pcolormesh(x, z, profile, shading="auto", cmap="magma")
    fig.colorbar(m, ax=ax, label="englacial concentration (kg/m^3)")
    ax.plot(x, Hc, "w", linewidth=0.5)
    ax.set_ylim(0, Hc.max() * 1.2 + 1)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("height above the bed (m)")

    ax = axes[2, 0]
    ax.plot(t, M_eng, label="englacial")
    ax.plot(t, M_sup, label="supraglacial")
    ax.plot(t, M_tot, "k", label="total")
    ax.set_ylabel("debris mass (kg)")
    ax.set_xlabel("time (years)")
    ax.legend()

    ax = axes[2, 1]
    ax.plot(t, F_in, label="input")
    ax.plot(t, F_melt, label="melt-out")
    ax.plot(t, F_out, label="output")
    ax.set_ylabel("mass flux (kg/yr)")
    ax.set_xlabel("time (years)")
    ax.legend()
    ax.set_title(f"conservation error at the end: {float(err[-1]):.3g} kg")

    fig.tight_layout()
    fig.savefig(figure, dpi=120)


if __name__ == "__main__":
    main(*sys.argv[1:4])
