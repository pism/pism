# Idealized debris-covered valley glacier

This example reproduces the setup of

> Y. Verhaegen and P. Huybrechts (2026). Coupling Debris Transport to 3D Higher-Order Ice
> Flow Dynamics to Model the Behavior and Climate Change Response of Debris-Covered
> Glaciers. *JGR Earth Surface*, 131, e2025JF008748. https://doi.org/10.1029/2025JF008748

(Sections 2.2 and 2.4): a 7.5 km long valley with a headwall merging into a linear bed
(slope 0.2) and a U-shaped cross section (100 m deep, 1 km wide), a linear mass balance
profile (`-0.003 (x - 1750 m)` m w.e./yr) and a debris source of `2e6` kg/yr on a
50 m x 350 m patch near the headwall, released 5 years into the run.

PISM's Blatter (higher-order) stress balance with pseudo-plastic sliding on a Mohr-Coulomb
yield stress and an isothermal Glen flow law computes the ice flow; the `transport` debris
model (`-debris transport`) buries the debris in the accumulation zone, advects it
englacially, melts it out in the ablation zone, moves it along the surface and downslope,
and removes it at the glacier margin. Its state (`debris_thickness`,
`englacial_debris_concentration`), the melt-out, input and removal rates, the debris-covered
fraction and the mass budget (scalar time series) are written to the output files.

Note that, in this version, the debris does not yet affect the surface mass balance
(the melt enhancement factor `ice_melt_enhancement` is computed and written but not applied),
so the glacier is the paper's *clean-ice* reference glacier and the debris a passive tracer.

## Ice dynamics and spin-up

PISM's Blatter solver cannot handle the thin, freshly nucleated ice of the first years at
grid spacings below about 100 m (its linear solve diverges, with or without the multigrid
preconditioner). `run.sh` therefore grows the glacier with the hybrid SSA+SIA stress
balance, using the same sliding law and basal strength, for `SPINUP` years (default 30)
and restarts the Blatter run from that state (warm-started with the SSA velocities). The
debris model runs during the spin-up as well, so the debris state and budget carry over.
Set `SPINUP=0` to run Blatter from the start (works at 100 m).

## Running

```
make                     # 25 m grid, 300 years, 4 processes: takes a while
make DX=100 DURATION=50  # a quick look

Note that `run.sh` only regenerates `input.nc`/`debris_input.nc` when they are missing:
run `make clean` (or delete them) before changing `DX`.
```

or step by step:

```
python3 create_input.py --dx 25 --debris-file debris_input.nc input.nc
N=4 DX=25 SPINUP=30 ./run.sh input.nc 300 debris
python3 plot_results.py debris_spatial.nc debris_scalar.nc figure.png
```

`run.sh` documents the full set of options. The marginal length scale of the debris removal
(`debris.transport.marginal_length_scale`) is set to the grid spacing, as in the paper.
