#!/usr/bin/env python

import PISM
import numpy as np
import pylab as plt
import ch_warming
year_length = float(365 * 86400)

z, times, T_ice, W_ice, T_ch, W_ch = ch_warming.run(Lz=200, Mz=201, T_final_years=4, R=20)

t = 86400 * np.linspace(0, 365, 366)
Ts = ch_warming.T_surface(t, mean=268.15, amplitude=6, summer_peak_day=365/2)
melt_season_length = np.count_nonzero(Ts > ch_warming.T_melting)
print("The melt season is {} weeks long.".format(melt_season_length / 7))

plt.figure(figsize=(10, 5))
plt.plot(t / 86400, Ts)
plt.xlabel("time, days")
plt.ylabel("Air temperature, Kelvin")
plt.grid()
plt.savefig("air-temperature.png")

N = z.size - 16
plt.figure(figsize=(10, 5))
plt.title("Ice temperature at depth over time")
plt.pcolormesh(times / year_length, z[N:], T_ice.T[N:, :])
plt.xlabel("time, years")
plt.ylabel("z, meters")
plt.colorbar()
plt.grid()
plt.savefig("ice-temperature.png")

N = z.size - 16
plt.figure(figsize=(10, 5))
plt.title("Temperature in the CH system")
plt.pcolormesh(times / year_length, z[N:], T_ch.T[N:, :])
plt.xlabel("time, years")
plt.ylabel("z, meters")
plt.colorbar()
plt.grid()
plt.savefig("ch-temperature.png")

N = z.size - 16
plt.figure(figsize=(10, 5))
plt.title("Water fraction in the CH system")
plt.pcolormesh(times / year_length, z[N:], W_ch.T[N:, :])
plt.xlabel("time, years")
plt.ylabel("z, meters")
plt.colorbar()
plt.grid()
plt.savefig("ch-water-fraction.png")

plt.figure(figsize=(10, 5))
N = z.size-16
for k in np.r_[N:z.size]:
    plt.plot(times/year_length, T_ice[:, k])

plt.xlabel("time, years")
plt.ylabel("temperature, Kelvin")
plt.title("Ice temperature at different depths")
plt.grid()
plt.savefig("ice-temperature-curves.png")
