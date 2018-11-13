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

            v[i, j].u = -y * radial_velocity
            v[i, j].v = x * radial_velocity

    v.update_ghosts()


def quiver(v, **kwargs):
    a = v.numpy()
    plt.quiver(a[:, :, 0], a[:, :, 1], **kwargs)


class MassTransport(object):

    def __init__(self, grid, part_grid=False):

        self.grid = grid

        if part_grid:
            grid.ctx().config().set_boolean("geometry.part_grid.enabled", True)

        self.v = PISM.model.create2dVelocityVec(grid)

        self.Q = PISM.IceModelVec2Stag()
        self.Q.create(grid, "Q", PISM.WITHOUT_GHOSTS)

        self.v_bc_mask = PISM.IceModelVec2Int()
        self.v_bc_mask.create(grid, "v_bc_mask", PISM.WITHOUT_GHOSTS)

        self.H_bc_mask = PISM.IceModelVec2Int()
        self.H_bc_mask.create(grid, "H_bc_mask", PISM.WITHOUT_GHOSTS)

        self.ge = PISM.GeometryEvolution(grid)

        geometry = PISM.Geometry(grid)
        self.geometry = geometry

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


def test():
    ctx = PISM.Context()
    Mx = 50
    My = 50
    grid = PISM.IceGrid_Shallow(ctx.ctx, 1, 1, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    import pylab as plt
    plt.rcParams['figure.figsize'] = (8.0, 8.0)

    mt = MassTransport(grid)

    mt.reset()
    levels = np.linspace(0, 1, 11)

    plt.subplot(2, 2, 1)
    mt.plot_thickness(levels, "time=0")

    mt.step(0.333)

    plt.subplot(2, 2, 2)
    mt.plot_thickness(levels, "time=0.333")

    mt.step(0.333)

    plt.subplot(2, 2, 3)
    mt.plot_thickness(levels, "time=0.666")

    mt.step(0.333)

    plt.subplot(2, 2, 4)
    mt.plot_thickness(levels, "time=0.999")

    plt.show()


if __name__ == "__main__":
    test()
