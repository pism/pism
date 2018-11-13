import PISM
from PISM.util import convert
import numpy as np
import pylab as plt

ctx = PISM.Context()
unit_system = ctx.unit_system
config = ctx.config

config.set_string("grid.ice_vertical_spacing", "equal")

k = config.get_double("constants.ice.thermal_conductivity")
T_melting = config.get_double("constants.fresh_water.melting_point_temperature")

EC = ctx.enthalpy_converter
pressure = np.vectorize(EC.pressure)
enthalpy = np.vectorize(EC.enthalpy)
temperature = np.vectorize(EC.temperature)
melting_temperature = np.vectorize(EC.melting_temperature)
water_fraction = np.vectorize(EC.water_fraction)


class EnthalpyColumn(object):
    "Set up the grid and arrays needed to run column solvers"

    def __init__(self, prefix, Mz, dt, Lz=1000.0):
        self.Lz = Lz
        self.z = np.linspace(0, self.Lz, Mz)

        param = PISM.GridParameters()
        param.Lx = 1e5
        param.Ly = 1e5
        param.z = PISM.DoubleVector(self.z)
        param.Mx = 3
        param.My = 3
        param.Mz = Mz

        param.ownership_ranges_from_options(1)

        self.dt = dt

        self.grid = PISM.IceGrid(ctx.ctx, param)
        grid = self.grid

        self.enthalpy = PISM.model.createEnthalpyVec(grid)

        self.strain_heating = PISM.model.createStrainHeatingVec(grid)
        self.strain_heating.set(0.0)

        self.u, self.v, self.w = PISM.model.create3DVelocityVecs(grid)
        self.u.set(0.0)
        self.v.set(0.0)
        self.w.set(0.0)

        self.sys = PISM.enthSystemCtx(grid.z(), prefix,
                                      grid.dx(), grid.dy(), self.dt,
                                      config,
                                      self.enthalpy,
                                      self.u, self.v, self.w,
                                      self.strain_heating,
                                      EC)

    def init(self):
        ice_thickness = self.Lz
        self.sys.init(1, 1, False, ice_thickness)


def ch_heat_flux(T_ice, T_ch, k, R):
    "Heat flux from the CH system into the ice. See equation 1."
    return (k / R**2) * (T_ch - T_ice)


def T_surface(time, mean, amplitude, summer_peak_day):
    "Surface temperature (cosine yearly cycle)"
    day_length = 86400
    summer_peak = summer_peak_day * day_length
    year_length = float(365 * day_length)

    t = np.mod(time - summer_peak, year_length) / year_length

    return mean + amplitude * np.cos(2 * np.pi * t)


def T_steady_state(H, z, T_surface, G, ice_k):
    "Steady state ice temperature given T_surface and the geothermal flux G."
    return T_surface + (H - z) * (G / ice_k)


def E_steady_state(H, z, T_surface, G, ice_k):
    return enthalpy(T_steady_state(H, z, T_surface, G, ice_k),
                    0.0, pressure(H - z))


def run(T_final_years=10.0, dt_days=1, Lz=1000, Mz=101, R=20, omega=0.005):
    """Run the one-column cryohydrologic warming setup"""

    T_final = convert(T_final_years, "years", "seconds")
    dt = convert(dt_days, "days", "seconds")

    ice = EnthalpyColumn("energy.enthalpy", Mz, dt, Lz=Lz)
    ch = EnthalpyColumn("energy.ch_warming", Mz, dt, Lz=Lz)

    H = ice.Lz

    z = np.array(ice.sys.z())
    P = pressure(H - z)

    z_coarse = np.array(ice.grid.z())
    P_coarse = pressure(H - z_coarse)

    T_mean_annual = 268.15    # mean annual temperature, Kelvin
    T_amplitude = 6         # surface temperature aplitude, Kelvin
    summer_peak_day = 365/2
    G = 0.0       # geothermal flux, W/m^2

    E_initial = E_steady_state(H, z_coarse, T_mean_annual, G, k)

    with PISM.vec.Access(nocomm=[ice.enthalpy, ice.strain_heating, ice.u, ice.v, ice.w,
                                 ch.enthalpy, ch.strain_heating, ch.u, ch.v, ch.w]):

        # set initial state for the ice enthalpy
        ice.enthalpy.set_column(1, 1, E_initial)

        # set initial state for the CH system enthalpy
        ch.enthalpy.set_column(1, 1, E_initial)

        times = []
        T_ice = []
        T_ch = []
        W_ice = []
        W_ch = []
        t = 0.0
        while t < T_final:
            times.append(t)

            # compute the heat flux from the CH system into the ice
            T_s = T_surface(t, T_mean_annual, T_amplitude, summer_peak_day)

            T_column_ice = temperature(ice.enthalpy.get_column_vector(1, 1), P_coarse)
            T_column_ch = temperature(ch.enthalpy.get_column_vector(1, 1), P_coarse)

            if R > 0:
                Q = ch_heat_flux(T_column_ice, T_column_ch, k, R)
            else:
                Q = np.zeros_like(P_coarse)

            # add it to the strain heating term
            ice.strain_heating.set_column(1, 1, Q)
            ch.strain_heating.set_column(1, 1, -Q)

            if T_s > T_melting:
                E_s = EC.enthalpy(T_melting, 0.0, 0.0)
            else:
                E_s = EC.enthalpy(T_s, 0.0, 0.0)

            # set boundary conditions and update the ice column
            ice.init()
            ice.sys.set_surface_dirichlet_bc(E_s)
            ice.sys.set_basal_heat_flux(G)
            x = ice.sys.solve()
            T_ice.append(temperature(x, P))
            W_ice.append(water_fraction(x, P))
            ice.sys.fine_to_coarse(x, 1, 1, ice.enthalpy)

            if T_s > T_melting:
                E_s = EC.enthalpy(T_melting, omega, 0.0)
                # Re-set enthalpy in the CH system during the melt season
                ch.enthalpy.set_column(1, 1, enthalpy(melting_temperature(P), omega, P))
            else:
                E_s = EC.enthalpy(T_s, 0.0, 0.0)

            # set boundary conditions and update the CH column
            ch.init()
            ch.sys.set_surface_dirichlet_bc(E_s)
            ch.sys.set_basal_heat_flux(G)
            x = ch.sys.solve()
            T_ch.append(temperature(x, P))
            W_ch.append(water_fraction(x, P))
            ch.sys.fine_to_coarse(x, 1, 1, ch.enthalpy)

            t += dt

    return z, np.array(times), np.array(T_ice), np.array(W_ice), np.array(T_ch), np.array(W_ch)
