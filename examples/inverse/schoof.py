class SchoofSSA1dExact:

    """
    schoof_ssa_exact
    Returns certain exact solutions of the ssa equation found in
    Schoof, A Variational Approach to Ice Stream Flow, 2006, pp 237-238.

    The PDE is:

    -d/dy (Bh/2 |1/2 du/dy|^(-2/3) du/dy) = f (1-|y/L|^m)

    on the domain -3L <= y <= 3L with periodic boundary conditions.  The resulting $u$
    is the downstream velocity on an infinite slab.

    """

    def __init__(self, L, m, B=1, h=1, f=1):
        self.L = float(L)
        self.m = float(m)
        self.f = float(f)
        self.h = float(h)
        self.B = float(B)
        self.scale = 2 * (f / (B * h)) ** 3.0

    def eval(self, x):
        L = self.L
        m = self.m
        W = (m + 1.0) ** (1.0 / m)
        u = abs(x / self.L)
        if u > W:
            v = 0
        else:
            v = -(L**4) * (
                (u**4 - (m + 1) ** (4.0 / m)) / 4
                - 3 * (u ** (m + 4) - (m + 1) ** (1 + 4.0 / m)) / ((m + 1) * (m + 4))
                + 3
                * (u ** (2 * m + 4) - (m + 1) ** (2 + 4.0 / m))
                / ((m + 1) ** 2 * (2 * m + 4.0))
                - (u ** (3 * m + 4) - (m + 1) ** (3 + 4.0 / m))
                / ((m + 1) ** 3 * (3 * m + 4))
            )

        v *= self.scale
        return v
