# dEBM-enhanced: implementation notes

This documents the **PISM implementation** of the `debm_enhanced` surface model — in
particular the internally-computed terrain-shaded insolation. For the underlying math
(ported from [`solshade`](https://github.com/amanchokshi/solshade)) see
[`solstis.md`](./solstis.md); this file describes what the C++ code actually does.

## Overview

`debm_enhanced` is dEBM-simple with the insolation that drives the insolation-melt term
replaced by a terrain-shaded surface-insolation field computed internally from the ice
surface elevation. At initialization it builds the terrain horizon, surface normals, and
sky-view factor; each time step it computes the daily terrain-shaded insolation.

Everything else (temperature/offset melt, albedo, refreezing, snow bookkeeping, the
atmospheric transmissivity, the empirical melt factors) is inherited unchanged from
dEBM-simple (`surface.debm_simple.*`).

## Source files

| File | Role |
|---|---|
| `src/coupler/surface/DEBMEnhanced.{hh,cc}` | the surface model; owns the engine and drives it from the dEBM-simple seams |
| `src/coupler/surface/TerrainInsolation.{hh,cc}` | engine: gathers the DEM, builds horizon / normals / sky-view factor, computes the daily insolation |
| `src/coupler/surface/terrain_insolation_kernel.{hh,cc}` | pure (PISM-free) numerical kernels: bilinear sampling, ray-march horizon, surface normal, sun position, sky-view factor |

`DEBMEnhanced` inherits from `DEBMSimple`. The two dEBM-simple override points it uses are
`update_insolation_input()` (once per time step) and `insolation_energy_series()` (per cell).

## Conventions

- Grid indexing matches PISM: `i` is the x-index (increasing **east**), `j` is the y-index
  (increasing **north**); the gathered global DEM is row-major `dem[j*Mx + i]`.
- **Azimuth is clockwise from north** (0 = N, π/2 = E), so a horizontal direction is
  `(East, North) = (sin A, cos A)`. The same convention is used for the horizon ray-march and
  the Sun azimuth, which makes the shadow test self-consistent. (The offline solshade file
  stores the azimuth axis as CCW-from-east; the two relate by `θ_file = 90° − A`.)
- Horizon angles are stored in **radians** internally; the `horizon` diagnostic is in radians.
- Latitude (`geometry.latitude`, degrees north) is converted to radians where needed.

## Pipeline

### One-time setup — `TerrainInsolation::init()` (called from `DEBMEnhanced::init_impl`)

1. **Gather the DEM to every rank.** The horizon ray-march reaches several km from each cell,
   crossing MPI subdomain boundaries. The DEM is gathered onto rank 0 (`put_on_proc0`) and
   broadcast (`MPI_Bcast`) so every rank holds the full global surface elevation
   (`dem[j*Mx+i]`, ~`Mx·My·8` bytes). Each rank then ray-marches only its **owned** cells
   against this shared array — embarrassingly parallel, no ghost communication, correct at
   domain edges.
2. **Surface normals** (per owned cell): upward unit normal from centred differences of the
   global DEM (one-sided at the domain edge): `nE = −∂z/∂E`, `nN = −∂z/∂N`, `nU = 1`,
   normalized (`surface_normal`).
3. **Horizon map** `horizon(azimuth, y, x)` (an `array::Array3D` with `n_directions` levels):
   for each azimuth `A_k = 2πk/n_directions`, `ray_horizon` marches in `step` increments to
   `max_distance`, bilinearly sampling the DEM, and takes the maximum elevation angle
   `atan2(z(d) − z₀, d)` along the ray. Rays stop at the domain boundary.
4. **Sky-view factor** `sky_view_factor(y, x)` (if `use_sky_view_factor`): the slope-corrected
   Dozier & Frew (1990) integral over the horizon column, using the cell's slope/aspect
   (recovered from the normal). Stored as an `array::Scalar` in `[0, 1]`.

### Per–time-step — `DEBMEnhanced::update_insolation_input(t, dt, …)`

1. **Refresh the horizon** if at least `update_interval` has elapsed since the last horizon
   computation (the surface elevation evolves): re-run `TerrainInsolation::init` on the
   current `geometry.ice_surface_elevation`. `update_interval = 0` recomputes every step.
2. **Solar parameters for the representative day** at `t + dt/2`, using PISM's analytic
   present-day orbit (`DEBMSimplePointwise::solar_declination_present_day` and
   `distance_factor_present_day`; `distance_factor = (d̄/d)² = 1/R²`). No ephemeris/Skyfield.
3. **Daily insolation field** `TerrainInsolation::daily_insolation(declination,
   distance_factor, latitude, result)` (see below), stored as the daily-mean `insolation`
   rate (`W m⁻²`, matching dEBM-simple's `insolation` diagnostic units).

### Per cell — `DEBMEnhanced::insolation_energy_series(i, j, …)`

The daily-mean rate is converted to the energy reaching the surface during each dEBM-simple
sub-step of length `dt_sub` (seconds): `result[k] = insolation(i, j) · dt_sub`. The downstream
dEBM-simple melt core (`melt_from_insolation`, which applies `atmosphere_transmissivity` and
the melt factors) consumes this exactly as it consumes dEBM-simple's analytic insolation.

## The daily insolation integral

`daily_insolation` integrates the diurnal cycle for the day's fixed declination `δ` and
distance factor `df`, then divides by the length of the day to a **daily-mean rate**
(`W m⁻²`). Because the integral runs over a full day (hour angle `H` swept over `[−π, π)` with
`M = 86400 / ephemeris_dt` midpoint samples), it is **independent of longitude**. For each
sample the Sun altitude `a` and azimuth `A` come from `sun_position`
(`sin a = sin φ sin δ + cos φ cos δ cos H`, standard azimuth), and the surface flux is split
into direct and diffuse:

```
I(i,j) = (1/T) · Σ_samples [ (1 − f) · S₀·df · max(0, n·s) · 1[a>0 and a>horizon(A)]   (direct)
                           +     f  · S₀·df · sin a · SVF · 1[a>0]            ] · Δt    (diffuse)
```

where `T` = 86400 s (so the leading `(1/T)·Σ…Δt` is the daily average).

- `S₀` = `surface.debm_simple.solar_constant` (1361 W m⁻²); `s` = ENU Sun unit vector;
  `n·s` = cosine of incidence on the tilted surface.
- `f` = `surface.debm_enhanced.diffuse_fraction` (the diffuse share); `SVF` = the cell's
  sky-view factor. The **diffuse term reaches shadowed cells** (it is gated only by `a > 0`,
  not the horizon), while the **direct term requires the Sun above the local horizon**.
- When `use_sky_view_factor` is off, `f` is forced to 0 and the result is **pure direct
  beam** — identical to the computation before the sky-view factor was added. `f = 0` gives
  the same.

This is a **top-of-atmosphere** decomposition: atmospheric attenuation is applied later by the
dEBM-simple melt core, not here. `f` is a fixed fraction, not derived from cloud cover.

## Configuration parameters (`surface.debm_enhanced.*`)

| Parameter | Default | Meaning |
|---|---|---|
| `horizon.n_directions` | `360` | azimuth directions for the horizon ray-march |
| `horizon.max_distance` | `5000` m | maximum ray-march distance |
| `horizon.step` | `100` m | ray-march step length |
| `horizon.ephemeris_dt` | `600` s | sub-daily step for the diurnal insolation integral |
| `update_interval` | `10` (365-day) years | how often the horizon is recomputed from the evolving geometry (0 = every step) |
| `use_sky_view_factor` | `yes` | compute the sky-view factor and apply the diffuse split (off → pure direct beam) |
| `diffuse_fraction` | `0.2` | isotropic diffuse share `f`, scaled by the sky-view factor |

## Diagnostics

- `insolation` — the daily-mean terrain-shaded insolation rate driving the melt (`W m⁻²`).
  This replaces dEBM-simple's analytic `insolation` diagnostic (also `W m⁻²`, but a mean over
  the melt period rather than a daily mean).
- `horizon` — the terrain horizon map `(azimuth, y, x)`, radians.
- `sky_view_factor` — the sky-view factor `(y, x)`, in `[0, 1]` (only if `use_sky_view_factor`).

## Verification

- The C++ `ray_horizon` kernel reproduces the solshade reference horizon
  (`insolation_RGI2000-v7.0-C-01-04374.nc`) to within **3×10⁻⁴ degrees** over an interior
  crop (see `verify_horizon.py` / `verify_horizon_main.cc`), which also confirms the
  azimuth-convention relation `θ_file = 90° − A`.
- `sky_view_factor` passes analytic checks: flat unobstructed → 1; flat with uniform horizon
  `h` → `cos²h`; fully enclosed → 0.
- With `diffuse_fraction = 0` (or `use_sky_view_factor = no`) the daily integral reduces
  exactly to the validated direct-beam computation.

## Approximations and limitations

- **Top-of-atmosphere, no explicit atmosphere here.** Attenuation is the dEBM-simple
  transmissivity applied downstream; there is no cloud/aerosol model. `diffuse_fraction` is a
  fixed constant, not a clearness-derived quantity.
- **Analytic present-day Sun position**, not a precise ephemeris (solshade used Skyfield /
  JPL DE440). Expect a few-percent difference in daily insolation, larger near shadow
  boundaries and sunrise/sunset.
- **No terrain-reflected component.** Only direct + isotropic diffuse.
- **Representative-day declination**: a single declination/distance (interval midpoint) is
  used across each update; accurate for the sub-monthly steps used in practice.
- **Static-within-interval horizon**: the horizon is recomputed only every `update_interval`,
  so it lags the evolving geometry between recomputations.
- **Domain edges**: rays are clamped at the domain boundary (no NaN), so near-edge cells
  differ slightly from the NaN-producing offline pipeline.

## Credits

The horizon, surface-normal, and direct-beam flux algorithms are re-implementations of
[`solshade`](https://github.com/amanchokshi/solshade) by Aman Chokshi (MIT License, © 2025;
Chokshi et al., JOSS, [doi:10.21105/joss.09944](https://doi.org/10.21105/joss.09944)) — see
the per-function credits in `terrain_insolation_kernel.cc`. The sky-view factor follows
Dozier & Frew (1990), [doi:10.1109/36.58986](https://doi.org/10.1109/36.58986). The solar
declination / distance factor are PISM's own dEBM-simple analytic orbit (Liou, 2002). Full
references are in [`solstis.md`](./solstis.md).
