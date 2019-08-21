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
    a = v.numpy()
    plt.quiver(a[:, :, 0], a[:, :, 1], **kwargs)


class MassTransport(object):

    def __init__(self, grid, part_grid=False):

        self.grid = grid

        if part_grid:
            grid.ctx().config().set_boolean("geometry.part_grid.enabled", True)

        self.v = PISM.IceModelVec2V(grid, "velocity", PISM.WITHOUT_GHOSTS)
        self.Q = PISM.IceModelVec2Stag(grid, "Q", PISM.WITHOUT_GHOSTS)
        self.v_bc_mask = PISM.IceModelVec2Int(grid, "v_bc_mask", PISM.WITHOUT_GHOSTS)
        self.H_bc_mask = PISM.IceModelVec2Int(grid, "H_bc_mask", PISM.WITHOUT_GHOSTS)

        self.ge = PISM.GeometryEvolution(grid)

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

    def plot_thickness(self, levels, title):
        import pylab as plt
        cm = plt.contour(self.grid.x(), self.grid.y(),
                         self.geometry.ice_thickness.numpy(), levels=levels)
        plt.clabel(cm)
        plt.grid()
        plt.title(title)

    def step(self, t_final, C=1):
        geometry = self.geometry
        t = 0.0
        j = 0
        while t < t_final:

            dt = PISM.max_timestep_cfl_2d(geometry.ice_thickness,
                                          geometry.cell_type,
                                          self.v).dt_max.value() * C

            if t + dt > t_final:
                dt = t_final - t

            log.message(2, "{}, {}\n".format(t, dt))

            self.ge.flow_step(geometry, dt,
                              self.v,
                              self.Q,
                              self.v_bc_mask,
                              self.H_bc_mask)

            geometry.ice_thickness.add(1.0, self.ge.thickness_change_due_to_flow())
            geometry.ice_area_specific_volume.add(1.0, self.ge.area_specific_volume_change_due_to_flow())
            geometry.ensure_consistency(0.0)

            t += dt
            j += 1

def volume(geometry):
    cell_area = geometry.ice_thickness.grid().cell_area()
    volume = ((geometry.ice_thickness.numpy() + geometry.ice_area_specific_volume.numpy()) *
              cell_area)
    return volume.sum()

def test():
    ctx = PISM.Context()
    Mx = int(ctx.config.get_double("grid.Mx"))
    My = int(ctx.config.get_double("grid.My"))
    grid = PISM.IceGrid_Shallow(ctx.ctx, 1, 1, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    import pylab as plt
    plt.rcParams['figure.figsize'] = (8.0, 8.0)

    mt = MassTransport(grid)

    mt.reset()
    levels = np.linspace(0, 1, 11)

    t = 0
    dt = 0.333
    C = 1.0

    for j in [1, 2, 4, 3]:
        plt.subplot(2, 2, j)
        mt.plot_thickness(levels, "time={}, V={}".format(t, volume(mt.geometry)))

        mt.step(dt, C)
        t += dt

    plt.show()


if __name__ == "__main__":
    test()
