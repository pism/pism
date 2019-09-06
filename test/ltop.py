#!/usr/bin/env python

import numpy as np

class LTOP(object):
    "Linear Theory of Orographic Precipitation (LTOP) model"

    tau_c = 1000.0
    "conversion time [s]"

    tau_f = 1000.0
    "fallout time [s]"

    P0 = 0.0
    'Background precipitation rate [mm hr-1]'

    P_scale = 1.0
    'Precipitation scale factor'

    Nm = 0.005
    'moist stability frequency [s-1]'

    Hw = 2500
    'Water vapor scale height [m]'

    latitude = 45.0
    "Latitude used to compute the Coriolis force"

    direction = 270.0
    "Wind direction, 0 is north, 270 is west"

    speed = 15.0
    "Wind speed"

    f = None
    "Coriolis force"

    u = None
    "u component of the wind velocity"

    v = None
    "v component of the wind velocity"

    Cw = None
    "Uplift sensitivity factor [kg m-3]"

    Theta_m = -6.5
    "moist adiabatic lapse rate [K / km]"

    rho_Sref = 7.4e-3
    "reference density [kg m-3]"

    gamma = -5.8
    "adiabatic lapse rate [K / km]"

    def __init__(self):
        self.update()

    def run(self, orography, dx, dy, truncate=True):
        "Compute orographic precipitation in mm/hour."
        # make sure derived constants are up to date
        self.update()

        eps = 1e-18

        nrows, ncols = orography.shape

        pad = max(nrows, ncols)

        h = np.pad(orography, pad, 'constant')
        nrows, ncols = h.shape

        h_hat = np.fft.fft2(h)

        x_freq = np.fft.fftfreq(ncols, dx / (2 * np.pi))
        y_freq = np.fft.fftfreq(nrows, dy / (2 * np.pi))

        kx, ky = np.meshgrid(x_freq, y_freq)

        # Intrinsic frequency sigma = U*k + V*l
        u0 = self.u
        v0 = self.v

        # $\sigma = U k + V l$, see paragraph after eq 5.
        sigma = u0 * kx + v0 * ky

        denominator = sigma**2 - self.f**2
        denominator[np.logical_and(np.fabs(denominator) < eps, denominator >= 0)] = eps
        denominator[np.logical_and(np.fabs(denominator) < eps, denominator  < 0)] = -eps

        m_squared = (self.Nm**2 - sigma**2) * (kx**2 + ky**2) / denominator

        m = np.sqrt(np.array(m_squared, dtype=np.complex))

        # Regularization
        nonzero = np.logical_and(m_squared >= 0, sigma != 0)
        m[nonzero] *= np.sign(sigma[nonzero])

        P_hat = h_hat * (self.Cw * 1j * sigma /
                         ((1 - 1j * m * self.Hw) *
                          (1 + 1j * sigma * self.tau_c) *
                          (1 + 1j * sigma * self.tau_f)))

        # Convert from wave domain back to space domain
        P = np.real(np.fft.ifft2(P_hat))

        # Remove padding
        if pad > 0:
            P = P[pad:-pad, pad:-pad]

        # convert to mm hr-1
        P *= 3600

        # Add background precipitation
        P += self.P0

        # Truncate
        if truncate:
            P[P < 0] = 0.0

        P *= self.P_scale

        return P

    def update(self):
        "Update derived constants"

        self.f = 2 * 7.2921e-5 * np.sin(self.latitude * np.pi / 180.0)

        self.u = -np.sin(self.direction * 2 * np.pi / 360) * self.speed
        self.v = -np.cos(self.direction * 2 * np.pi / 360) * self.speed

        self.Cw = self.rho_Sref * self.Theta_m / self.gamma

def triangle_ridge_grid(dx=5e4, dy=5e4):
    "Allocate the grid for the synthetic geometry test."
    x_min, x_max   = -100e3, 100e3
    y_min, y_max   = -100e3, 100e3

    Mx = int((x_max - x_min) / dx) + 1
    My = int((y_max - y_min) / dy) + 1

    x = np.linspace(x_min, x_max, Mx)
    y = np.linspace(y_min, y_max, My)

    return x, dx, y, dy

def triangle_ridge(x, A=500.0, d=50e3):
    "Create the 'triangle ridge' topography"
    return np.maximum(A * (1 - np.fabs(x) / d), 0)

def triangle_ridge_exact(x, u, Cw, tau, A=500.0, d=50e3):
    assert d > 0

    C = Cw * u * A / d
    Ut = u * tau

    xc = Ut * np.log(2 - np.exp(-d / Ut))

    def P(x):
        if x < 0 and x >= -d:
            return C * (1.0 - np.exp(-(x + d) / Ut))
        elif x >= 0 and x <= xc:
            return C * (np.exp(-x / Ut) * (2 - np.exp(-d / Ut)) - 1)
        else:
            return 0

    try:
        return 3600 * np.array([P(t) for t in x])
    except TypeError:
        return 3600 * P(x)

def max_error(spacing, direction):
    model = LTOP()
    # Set conversion time to zero (we could set fallout time to zero instead: it does not
    # matter which one is zero)
    model.tau_c = 0.0
    model.Hw = 0.0
    model.direction = direction
    model.latitude = 0.0

    if direction == 90 or direction == 270:
        # east or west
        x, dx, y, dy = triangle_ridge_grid(dx=spacing)
        t = x

        h = triangle_ridge(t)
        orography = np.tile(h, (len(y), 1))

        P = model.run(orography, dx, dy)
        P = P[len(y) // 2, :]
    else:
        # north or south
        x, dx, y, dy = triangle_ridge_grid(dy=spacing)
        t = y

        h = triangle_ridge(t)
        orography = np.tile(h, (len(x), 1)).T

        P = model.run(orography, dx, dy)
        P = P[:, len(x) // 2]

    if direction == 0 or direction == 90:
        P_exact = triangle_ridge_exact(-t, model.speed, model.Cw, model.tau_f)
    else:
        P_exact = triangle_ridge_exact(t,  model.speed, model.Cw, model.tau_f)

    return np.max(np.fabs(P - P_exact))

def convergence_rate(error, direction, plot):
    dxs = [2000, 1000]
    errors = [error(dx, direction) for dx in dxs]

    p = np.polyfit(np.log10(dxs), np.log10(errors), 1)

    if plot:
        import pylab as plt

        def log_plot(x, y, style, label):
            plt.plot(np.log10(x), np.log10(y), style, label=label)
            plt.xticks(np.log10(x), x)

        def log_fit_plot(x, p, label):
            plt.plot(np.log10(x), np.polyval(p, np.log10(x)), label=label)

        wind_direction = {0 : "north", 90 : "east", 180 : "south", 270 : "west"}

        plt.figure()
        plt.title("Precipitation errors (wind direction: {})".format(wind_direction[direction]))
        log_fit_plot(dxs, p, "polynomial fit (dx^{:1.4})".format(p[0]))
        log_plot(dxs, errors, 'o', "errors")
        plt.legend()
        plt.grid()
        plt.xlabel("grid spacing (meters)")
        plt.ylabel("log10(error)")
        plt.show()

    return p[0]

def ltop_test(plot=False):
    "Comparing to the 'triangle ridge' exact solution"
    assert convergence_rate(max_error,   0, plot) > 1.9
    assert convergence_rate(max_error,  90, plot) > 1.9
    assert convergence_rate(max_error, 180, plot) > 1.9
    assert convergence_rate(max_error, 270, plot) > 1.9

if __name__ == "__main__":
    ltop_test(plot=True)
