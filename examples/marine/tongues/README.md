### Testing eigen-calving in a setup with "tongues"

`generate_input.py` creates an input file with several "outlet
glaciers" of various widths, starting from 1, 2, ... grid points.

Run `make default` to run without calving, `make eigen_calving` to
turn on eigen-calving (seems to be unstable on its own because some
tongues develop very thin "offshoots" that are only weakly connected
to grounded ice).

Run `make eigen_plus_thickness_calving` to run with both eigen-calving
and thickness-calving on. This setup shows that one- and two-cell-wide
tongues *do not calve* at all in this setup.

This justifies the removal of one-cell-wide tongues using an
additional geometric (and strain-rate-independent) calving rule.

I'm not sure what to do with tongues that are two cells wide. I think
they will require some changes in the eigen-calving code.

#### Side note:

Strictly speaking, the eigen-calving model should produce a **0**
(zero) calving rate in this setup.

Recall that principal strain rates are computed as follows:

    const PetscScalar A = 0.5 * (u_x + v_y),  // A = (1/2) trace(D)
      B   = 0.5 * (u_x - v_y),
      Dxy = 0.5 * (v_x + u_y),  // B^2 = A^2 - u_x v_y
      q   = sqrt(PetscSqr(B) + PetscSqr(Dxy));
    eigen1 = A + q;
    eigen2 = A - q; // q >= 0 so e1 >= e2

In this setup there is no variation in ice thickness in each "tongue",
so `taud_x = 0`, `u_x = v_x = 0` and `u = u_y = 0`, so the code above
should be equivalent to

    eigen1 = v_y;
    eigen2 = 0; // q >= 0 so e1 >= e2

Now, the calving rule is

    if (eigen2 > eigenCalvOffset &&
        eigen1 > 0.0) { // if spreading in all directions
      calving_rate_horizontal = m_K * eigen1 * (eigen2 - eigenCalvOffset);
    } else {
      calving_rate_horizontal = 0.0;
    }

and `eigenCalvOffset == 0`, so `calving_rate_horizontal` should be zero.
