import numpy as np
import PISM

"""This is a half-baked test for the mass transport code that sets up
a rotational velocity field and then steps forward in time for one revolution.

It is not ready to be used as a regression or a verification test. All
I can see right now is that our mass transport scheme is very
diffusive, which is not a surprise."""

log = PISM.Context().log

def disc(thickness, x0, y0, H, R):
    """Set ice thickness to H within the disc centered at (x0,y0) of
    radius R and 0 elsewhere.
    """

    grid = thickness.grid()

    with PISM.vec.Access(nocomm=thickness):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            d2 = (x - x0)**2 + (y - y0)**2
            r2 = R**2
            if d2 <= r2:
                thickness[i, j] = H
            else:
                thickness[i, j] = 0.0

    thickness.update_ghosts()


def set_ice_thickness(output, time):
    """Exact solution at time time. Corresponds to a disc that is rotated
    around the origin (one revolution per time unit)."""

    grid = output.grid()

    L = min(grid.Lx(), grid.Ly())

    # 1 revolution per 1 time unit
    phi = 2 * np.pi * time

    M = np.matrix([[np.cos(phi), -np.sin(phi)],
                   [np.sin(phi),  np.cos(phi)]])

    # center coordinates at time 0
    x0, y0 = 0.5 * L, 0.0

    R = 0.25 * L
    x, y = M * np.matrix([x0, y0]).T

    disc(output, x, y, 1, R)


def set_velocity(v):
    """Initialize the velocity field to a rigid rotation around the
    origin.

    """
    grid = v.grid()

    radial_velocity = 2 * np.pi

    with PISM.vec.Access(nocomm=v):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)

            v[i, j].u =  y * radial_velocity
            v[i, j].v = -x * radial_velocity

    v.update_ghosts()


def quiver(v, **kwargs):
    a = v.to_numpy()
    plt.quiver(a[:, :, 0], a[:, :, 1], **kwargs)


class MassTransport(object):

    def __init__(self, grid, scheme="pism", part_grid=False):

        self.new_scheme = False
        self.grid = grid

        self.v = PISM.Vector(grid, "velocity")
        self.Q = PISM.Staggered(grid, "Q")
        self.H_bc_mask = PISM.Scalar(grid, "H_bc_mask")

        self.ge = PISM.GeometryEvolution(grid)

        if scheme == "pism":
            self.new_scheme = False

            if part_grid:
                grid.ctx().config().set_flag("geometry.part_grid.enabled", True)
        elif scheme == "mpdata":
            self.transport_scheme = PISM.MPDATA2(grid, 2)
            self.new_scheme = True
        else:
            schemes = {"upwind" : PISM.PISM_UNO_UPWIND1,
                       "lax-wendroff": PISM.PISM_UNO_LAX_WENDROFF,
                       "fromm": PISM.PISM_UNO_FROMM,
                       "uno2": PISM.PISM_UNO_2,
                       "uno3": PISM.PISM_UNO_3}

            self.transport_scheme = PISM.UNO(grid, schemes[scheme])
            self.new_scheme = True

        self.geometry = PISM.Geometry(grid)

        self.reset()

    def reset(self):
        geometry = self.geometry
        # grid info
        geometry.latitude.set(0.0)
        geometry.longitude.set(0.0)
        # environment
        geometry.bed_elevation.set(-10.0)
        geometry.sea_level_elevation.set(0.0)
        # ice
        set_ice_thickness(geometry.ice_thickness, time=1)
        geometry.ice_area_specific_volume.set(0.0)

        geometry.ensure_consistency(0.0)

        set_velocity(self.v)

        self.Q.set(0.0)

    def plot_thickness(self, ax, levels):

        cm = ax.contour(self.grid.x(), self.grid.y(),
                        self.geometry.ice_thickness.to_numpy(), levels=levels)
        ax.clabel(cm)
        ax.grid()

    def step(self, t_final, C=1):
        geometry = self.geometry
        t = 0.0
        j = 0
        while t < t_final:

            cfl = PISM.max_timestep_cfl_2d(geometry.ice_thickness,
                                           geometry.cell_type,
                                           self.v)
            dt = cfl.dt_max.value() * C

            if t + dt > t_final:
                dt = t_final - t

            log.message(3, "{}, {}\n".format(t, dt))

            if not self.new_scheme:
                self.ge.flow_step(geometry, dt,
                                  self.v,
                                  self.Q,
                                  self.H_bc_mask)

                geometry.ice_thickness.add(1.0, self.ge.thickness_change_due_to_flow())
                geometry.ice_area_specific_volume.add(1.0, self.ge.area_specific_volume_change_due_to_flow())
            else:
                self.transport_scheme.update(dt, geometry.cell_type, geometry.ice_thickness, self.v, True)
                geometry.ice_thickness.copy_from(self.transport_scheme.x())

            geometry.ensure_consistency(0.0)

            t += dt
            j += 1

def volume(geometry):
    cell_area = geometry.ice_thickness.grid().cell_area()
    return (PISM.sum(geometry.ice_thickness) + PISM.sum(geometry.ice_area_specific_volume)) * cell_area

def test():
    ctx = PISM.Context()
    Mx = int(ctx.config.get_number("grid.Mx"))
    My = int(ctx.config.get_number("grid.My"))
    grid = PISM.Grid.Shallow(ctx.ctx, 1, 1, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)

    opt = PISM.OptionString("-scheme", "mass transport scheme", "pism")

    mt = MassTransport(grid, opt.value())

    mt.reset()
    levels = np.linspace(0, 2, 41)[1:]

    t0 = 0
    tf = 1
    t = t0
    dt = 0.25
    C = PISM.OptionReal(ctx.unit_system, "-cfl", "CFL number", "1", 0.5).value()
    N = int(tf / dt)

    mt.plot_thickness(ax, levels)

    Range = mt.geometry.ice_thickness.range()
    print(f"time={t}, min={Range[0]}, max={Range[1]}, V={volume(mt.geometry)}")

    for j in range(N):
        mt.step(dt, C)
        t += dt

        mt.plot_thickness(ax, levels)

        Range = mt.geometry.ice_thickness.range()

        print(f"time={t}, min={Range[0]}, max={Range[1]}, V={volume(mt.geometry)}")

    plt.show()

    mt.geometry.ice_thickness.dump(ctx.config.get_string("output.file"))

if __name__ == "__main__":
    test()
