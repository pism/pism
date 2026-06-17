# Solar insolation math (`solstis`)

This note documents the mathematics behind the two computations we are porting to C/GPU
from [`solshade`](https://github.com/amanchokshi/solshade): the **horizon map**
(`compute_horizon_map`) and the **terrain-aware solar flux time series**
(`compute_flux_timeseries`). It describes only the math required for those two functions
and their immediate inputs.

All angles are in radians unless a quantity is explicitly labeled in degrees. The local
reference frame is **East–North–Up (ENU)**, and azimuth is measured **clockwise from true
north** (Skyfield's convention).

A couple things worth flagging for the C/GPU port:

- This is direct-beam only — no atmospheric attenuation, no diffuse/sky-view, no terrain-reflected component. If dEBM-enhanced needs those, they aren't in these two functions.
  - Azimuth follows Skyfield's clockwise-from-north convention, and the horizon map is stored in degrees while the flux kernel works in radians — easy to get wrong when porting.

---

## 1. Geometry inputs

Both computations operate on a digital elevation model (DEM) `z(y, x)` sampled on a regular
grid with spacing `dx`, `dy` (meters). Two derived geometric quantities are needed.

### 1.1 Surface slope, aspect, and the surface-normal vector

From the DEM gradients

```
dz/dy, dz/dx = gradient(z; dy, dx)
```

the slope and aspect angles are

$$
\text{slope} = \arctan\!\sqrt{\left(\tfrac{\partial z}{\partial x}\right)^2 + \left(\tfrac{\partial z}{\partial y}\right)^2},
\qquad
\text{aspect} = \operatorname{atan2}\!\left(-\tfrac{\partial z}{\partial x},\ \tfrac{\partial z}{\partial y}\right),
$$

with aspect wrapped into `[0, 2π)`. The **outward surface normal** in ENU components is

$$
\hat{\mathbf n} =
\begin{pmatrix} n_E \\ n_N \\ n_U \end{pmatrix}
=
\begin{pmatrix}
\sin(\text{slope})\,\sin(\text{aspect}) \\
\sin(\text{slope})\,\cos(\text{aspect}) \\
\cos(\text{slope})
\end{pmatrix},
$$

then normalized to unit length, $\hat{\mathbf n} \leftarrow \hat{\mathbf n}/\lVert\hat{\mathbf n}\rVert$.
A flat cell ($\text{slope}=0$) gives $\hat{\mathbf n} = (0,0,1)$, i.e. straight up.

### 1.2 Solar position as an ENU unit vector

The Sun's apparent topocentric **altitude** $a$, **azimuth** $A$, and **distance** $r$ (in AU)
are obtained from a high-precision ephemeris (Skyfield + JPL DE440), including light-time and
aberration but **without** atmospheric refraction:

```
apparent = (earth + observer).at(t).observe(sun).apparent()
r   = apparent.distance().au
a, A = apparent.altaz()
```

The unit vector pointing from the observer **toward the Sun** in ENU is

$$
\hat{\mathbf s} =
\begin{pmatrix} s_E \\ s_N \\ s_U \end{pmatrix}
=
\begin{pmatrix}
\cos a \,\sin A \\
\cos a \,\cos A \\
\sin a
\end{pmatrix}.
$$

Note $s_U = \sin a$ is positive only when the Sun is above the (astronomical) horizon.

---

## 2. `compute_horizon_map` — terrain horizon angles

**Goal.** For every DEM pixel, compute the elevation angle of the local terrain horizon as a
function of azimuth. This answers: *"In direction $\theta$, how high must the Sun be before it
clears the surrounding topography?"*

### 2.1 Ray casting

For each pixel $(y_0, x_0)$ with elevation $z_0$, and for each of $N_\text{dir}$ azimuthal
directions $\theta_k = 2\pi k / N_\text{dir}$, a ray is marched outward in steps of `step`
meters up to `max_distance`. The horizontal sample offsets (in pixel units) are

$$
\Delta x_\text{pix} = \frac{\cos\theta_k \; d}{dx}, \qquad
\Delta y_\text{pix} = \frac{\sin\theta_k \; d}{dy},
$$

for sample distances $d \in \{d_1, d_2, \dots\}$. Elevations $z(d)$ along the ray are obtained
by **bilinear interpolation** of the DEM (`map_coordinates`).

### 2.2 Horizon angle

At each sample the elevation angle subtended by the terrain relative to the origin pixel is

$$
\phi(d) = \operatorname{atan2}\big(z(d) - z_0,\; d\big),
$$

and the horizon angle in direction $\theta_k$ is the maximum over the ray:

$$
H(y_0, x_0, \theta_k) = \max_{d} \; \phi(d)
= \max_{d} \; \operatorname{atan2}\big(z(d) - z_0,\; d\big).
$$

If no valid (non-NaN) samples exist along a ray, the horizon there is NaN. The result is a
3‑D array $H(\theta, y, x)$ stored in **degrees**.

This is the standard ray-traced horizon / "hill-shading horizon angle" used in topographic
shortwave-radiation models (e.g. GRASS GIS `r.sun`, SAGA, TopoCalc).

### 2.3 Azimuth interpolation

Because the flux step queries the horizon at the *Sun's* azimuth (which need not coincide with
a sampled $\theta_k$), $H(\cdot, y, x)$ is interpolated over azimuth with a **periodic cubic
spline** (falling back to linear when coverage is sparse or the spline produces unphysical
values outside $[-90°, 90°]$). Query azimuths are wrapped into $[0, 360°)$ before evaluation.

---

## 3. `compute_flux_timeseries` — terrain-aware insolation

**Goal.** For each time $t$ and pixel $(y, x)$, compute the instantaneous top-of-atmosphere
shortwave flux on the sloped, possibly shadowed surface.

The model combines three effects:

1. **Cast-shadow test** (does terrain block the Sun?),
2. **Lambertian cosine projection** onto the tilted surface, and
3. **Inverse-square distance scaling** of the solar "constant."

### 3.1 Shadow test

The pixel is illuminated only if the Sun's altitude exceeds the local horizon angle in the
Sun's azimuth direction $A(t)$:

$$
\text{illuminated}(t, y, x) \iff a(t) > H\big(A(t), y, x\big).
$$

Otherwise the flux is 0.

### 3.2 Surface projection (cosine of incidence)

For an illuminated pixel, the cosine of the angle between the surface normal and the Sun
direction is the dot product of the two ENU unit vectors,

$$
\cos\xi = \hat{\mathbf n}\cdot\hat{\mathbf s}
= n_E s_E + n_N s_N + n_U s_U,
$$

clamped to be non-negative (a self-shadowed back-facing slope receives no direct beam):

$$
\mu = \max\!\big(0,\ \hat{\mathbf n}\cdot\hat{\mathbf s}\big).
$$

### 3.3 Distance scaling and total flux

The beam irradiance on the surface is the reference top-of-atmosphere flux $S_0$ scaled by the
inverse square of the Earth–Sun distance $r(t)$ (in AU, so $r\approx 1$):

$$
F(t, y, x) =
\begin{cases}
\dfrac{S_0}{r(t)^2}\;\mu(t,y,x), & a(t) > H\big(A(t), y, x\big),\\[1.2em]
0, & \text{shadowed},\\[0.4em]
\text{NaN}, & H \text{ is NaN (outside domain)}.
\end{cases}
$$

The default reference flux is the solar constant $S_0 = 1361\ \mathrm{W\,m^{-2}}$, and $r(t)$ is
clipped to a small positive minimum to avoid division blow-up. The result is a flux time series
$F(t, y, x)$ in $\mathrm{W\,m^{-2}}$.

This is direct-beam, terrain-resolved insolation only: it accounts for cast shadows, slope/aspect
orientation, and Earth–Sun distance, but **not** atmospheric attenuation, diffuse/sky-view, or
reflected (terrain-bounce) components.

---

## 4. Summary of the pipeline

```
DEM z(y,x) ──> slope, aspect ──> surface normal  n̂(y,x)          [§1.1]
DEM z(y,x) ──> ray casting   ──> horizon map      H(θ,y,x)        [§2]

ephemeris  ──> a(t), A(t), r(t) ──> sun vector    ŝ(t)            [§1.2]

F(t,y,x) = (S₀ / r²) · max(0, n̂·ŝ) · 1[ a > H(A,y,x) ]            [§3]
```

---

## References

- **Skyfield** (solar ephemeris): B. Rhodes, *Skyfield: High precision research-grade
  positions for planets and Earth satellites*, Astrophysics Source Code Library, 2019.
  <https://rhodesmill.org/skyfield/>
- **JPL DE440 ephemerides**: R. S. Park, W. M. Folkner, J. G. Williams, D. H. Boggs,
  *The JPL Planetary and Lunar Ephemerides DE440 and DE441*, The Astronomical Journal, 161:105,
  2021. doi:[10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414)
- **NumPy** (gradients, surface normals): C. R. Harris et al., *Array programming with NumPy*,
  Nature 585, 357–362, 2020. doi:[10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
- **TopoCalc** (sky-view / terrain horizon factors): S. Havens, *TopoCalc: Sky view and terrain
  configuration factors*, 2021. <https://github.com/USDA-ARS-NWRC/topocalc>
- **GRASS GIS `r.sun`** (ray-traced horizon & solar radiation): M. Neteler, M. H. Bowman,
  M. Landa, M. Metz, *GRASS GIS: A multi-purpose open source GIS*, Environmental Modelling &
  Software 31, 124–130, 2012. doi:[10.1016/j.envsoft.2011.11.014](https://doi.org/10.1016/j.envsoft.2011.11.014)
- **SAGA** (terrain analysis): O. Conrad et al., *System for Automated Geoscientific Analyses
  (SAGA) v. 2.1.4*, Geoscientific Model Development 8, 1991–2007, 2015.
  doi:[10.5194/gmd-8-1991-2015](https://doi.org/10.5194/gmd-8-1991-2015)
- **pvlib python** (irradiance modeling, related tool): W. F. Holmgren, C. W. Hansen,
  M. A. Mikofski, *pvlib python: a python package for modeling solar energy systems*,
  Journal of Open Source Software 3(29), 884, 2018.
  doi:[10.21105/joss.00884](https://doi.org/10.21105/joss.00884)
- **solshade** (source of the algorithms documented here): A. Chokshi et al., Journal of Open
  Source Software. doi:[10.21105/joss.09944](https://doi.org/10.21105/joss.09944)
</content>
</invoke>
