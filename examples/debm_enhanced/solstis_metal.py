# Copyright (C) 2026 PISM Authors
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""macmetalpy (Apple Metal) backends for the solstis terrain ops.

All Metal Shading Language kernel sources and their Python launch wrappers live
here, separate from the C/ctypes backends and the numpy/xarray glue in
``solstis.py``. macmetalpy's ``RawKernel`` takes the MSL source verbatim (no
``_MSL_HEADER`` injection), so each kernel below is self-contained.

Public functions:
    metal_horizon(dem, dx_pix, dy_pix, distances, batch_mb=512, pre=None)
    metal_upload(dem, distances)                  # pre-upload for kernel-only timing
    metal_slope(dem, dx, dy)
    metal_aspect(dem, dx, dy)
    metal_normals(dem, dx, dy)
"""
from __future__ import annotations

import os

import numpy as np

# The C dylib and macmetalpy each link their own libomp; loading both in one
# process trips OpenMP's "multiple runtimes" abort. They don't actually share
# OpenMP state here, so allow the duplicate. Must be set before macmetalpy import.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")


# --------------------------------------------------------------------------- #
# Horizon map (macmetalpy RawKernel)
# --------------------------------------------------------------------------- #
_MSL = r"""
#include <metal_stdlib>
using namespace metal;

constant float INVALID = -1.0e30f;

// Bilinear sample replicating scipy map_coordinates(order=1, constant, cval=NaN),
// but returning the INVALID sentinel instead of NaN (fast-math safe).
inline float bilinear(device const float* dem, int H, int W, float fy, float fx) {
    float ffy = floor(fy), ffx = floor(fx);
    int y0 = (int)ffy, x0 = (int)ffx;
    float ty = fy - ffy, tx = fx - ffx;
    float w00 = (1.0f-ty)*(1.0f-tx), w01 = (1.0f-ty)*tx;
    float w10 = ty*(1.0f-tx),        w11 = ty*tx;
    float acc = 0.0f;
    if (w00 > 0.0f) { if (y0<0||y0>=H||x0<0||x0>=W) return INVALID;       acc += w00*dem[y0*W+x0]; }
    if (w01 > 0.0f) { if (y0<0||y0>=H||x0+1<0||x0+1>=W) return INVALID;   acc += w01*dem[y0*W+(x0+1)]; }
    if (w10 > 0.0f) { if (y0+1<0||y0+1>=H||x0<0||x0>=W) return INVALID;   acc += w10*dem[(y0+1)*W+x0]; }
    if (w11 > 0.0f) { if (y0+1<0||y0+1>=H||x0+1<0||x0+1>=W) return INVALID; acc += w11*dem[(y0+1)*W+(x0+1)]; }
    return acc;
}

kernel void horizon_map(device const float* dem       [[buffer(0)]],
                        device const float* dx_pix    [[buffer(1)]],
                        device const float* dy_pix    [[buffer(2)]],
                        device const float* distances [[buffer(3)]],
                        device float*       out       [[buffer(4)]],
                        device const float* params    [[buffer(5)]],
                        uint gid [[thread_position_in_grid]]) {
    int H    = (int)params[0];
    int W    = (int)params[1];
    int ndir = (int)params[2];
    int ns   = (int)params[3];
    uint npix = (uint)(H * W);
    if (gid >= npix) return;

    int iy = (int)(gid / (uint)W);
    int ix = (int)(gid % (uint)W);

    float elev0 = bilinear(dem, H, W, (float)iy, (float)ix);
    bool bad0 = (elev0 < -1.0e29f);

    for (int idir = 0; idir < ndir; ++idir) {
        float best = INVALID;
        bool found = false;
        if (!bad0) {
            int base = idir * ns;
            for (int k = 1; k < ns; ++k) {
                float s = bilinear(dem, H, W,
                                   (float)iy + dy_pix[base + k],
                                   (float)ix + dx_pix[base + k]);
                if (s < -1.0e29f) continue;          // invalid sample
                float ang = atan2(s - elev0, distances[k]);
                if (!found || ang > best) { best = ang; found = true; }
            }
        }
        out[(uint)idir * npix + gid] =
            found ? best * (180.0f / 3.14159265358979323846f) : NAN;
    }
}
"""

_KERNEL = None


def _metal_kernel():
    global _KERNEL
    if _KERNEL is None:
        import macmetalpy as xp
        _KERNEL = xp.RawKernel(_MSL, "horizon_map")
    return _KERNEL


def metal_horizon(dem, dx_pix, dy_pix, distances, batch_mb=512, pre=None):
    """Run the Metal kernel, tiling over directions so no single output buffer
    exceeds ``batch_mb`` MB (this device caps a Metal buffer near ~1.5-2 GB, so
    the full DEM's 2.34 GB output must be split). ``pre`` reuses an already
    uploaded ``(d_dem, d_dist)`` -- used by --metal-exclude-transfer."""
    import macmetalpy as xp

    H, W = dem.shape
    n_dir, ns = dx_pix.shape
    kern = _metal_kernel()

    plane_bytes = H * W * 4
    db = max(1, int(batch_mb) * 1024 * 1024 // plane_bytes)   # dirs per tile
    db = min(db, n_dir)

    if pre is None:
        d_dem = xp.asarray(np.ascontiguousarray(dem, np.float32))
        d_dist = xp.asarray(np.ascontiguousarray(distances, np.float32))
    else:
        d_dem, d_dist = pre

    out = np.empty((n_dir, H, W), dtype=np.float32)
    for b0 in range(0, n_dir, db):
        b1 = min(b0 + db, n_dir)
        nb = b1 - b0
        d_dx = xp.asarray(np.ascontiguousarray(dx_pix[b0:b1].reshape(-1), np.float32))
        d_dy = xp.asarray(np.ascontiguousarray(dy_pix[b0:b1].reshape(-1), np.float32))
        d_par = xp.asarray(np.array([H, W, nb, ns], np.float32))
        d_out = xp.empty((nb, H, W), dtype=np.float32)
        kern(H * W, [d_dem, d_dx, d_dy, d_dist, d_out, d_par])
        out[b0:b1] = d_out.get().reshape(nb, H, W)
    return out


def metal_upload(dem, distances):
    import macmetalpy as xp

    return (
        xp.asarray(np.ascontiguousarray(dem, np.float32)),
        xp.asarray(np.ascontiguousarray(distances, np.float32)),
    )


# --------------------------------------------------------------------------- #
# Slope / aspect / normals (solshade compute_slope_aspect_normals)
#
# Local-stencil ops: a gradient via numpy.gradient's edge_order=1 rule (central
# differences interior, first-order one-sided on the border), then the GIS
# slope/aspect/ENU-normal math. Each is a standalone per-pixel kernel; the
# params buffer carries [H, W, dx, dy].
# --------------------------------------------------------------------------- #
_GRAD_MSL = r"""
#include <metal_stdlib>
using namespace metal;

// numpy.gradient (edge_order=1) at pixel (iy, ix); dzdx is d/dx (cols).
inline void grad(device const float* dem, int H, int W, int iy, int ix,
                 float dx, float dy, thread float& gx, thread float& gy) {
    int r = iy * W;
    if (W < 2)            gx = 0.0f;
    else if (ix == 0)     gx = (dem[r+1] - dem[r]) / dx;
    else if (ix == W-1)   gx = (dem[r+W-1] - dem[r+W-2]) / dx;
    else                  gx = (dem[r+ix+1] - dem[r+ix-1]) / (2.0f*dx);

    if (H < 2)            gy = 0.0f;
    else if (iy == 0)     gy = (dem[W+ix] - dem[ix]) / dy;
    else if (iy == H-1)   gy = (dem[(H-1)*W+ix] - dem[(H-2)*W+ix]) / dy;
    else                  gy = (dem[(iy+1)*W+ix] - dem[(iy-1)*W+ix]) / (2.0f*dy);
}

// Aspect in radians; pin the flat case to 0 (Metal's fast-math atan2(0,0) is
// NaN, whereas numpy's atan2(0,0) == 0).
inline float aspect_rad(float gx, float gy) {
    return (gx == 0.0f && gy == 0.0f) ? 0.0f : atan2(-gx, gy);
}
"""

_SLOPE_SRC = _GRAD_MSL + r"""
kernel void slope_map(device const float* dem    [[buffer(0)]],
                      device float*       out    [[buffer(1)]],
                      device const float* params [[buffer(2)]],
                      uint gid [[thread_position_in_grid]]) {
    int H = (int)params[0], W = (int)params[1];
    float dx = params[2], dy = params[3];
    if (gid >= (uint)(H*W)) return;
    int iy = (int)(gid / (uint)W), ix = (int)(gid % (uint)W);
    float gx, gy; grad(dem, H, W, iy, ix, dx, dy, gx, gy);
    out[gid] = atan(sqrt(gx*gx + gy*gy)) * (180.0f / 3.14159265358979323846f);
}
"""

_ASPECT_SRC = _GRAD_MSL + r"""
kernel void aspect_map(device const float* dem    [[buffer(0)]],
                       device float*       out    [[buffer(1)]],
                       device const float* params [[buffer(2)]],
                       uint gid [[thread_position_in_grid]]) {
    int H = (int)params[0], W = (int)params[1];
    float dx = params[2], dy = params[3];
    if (gid >= (uint)(H*W)) return;
    int iy = (int)(gid / (uint)W), ix = (int)(gid % (uint)W);
    float gx, gy; grad(dem, H, W, iy, ix, dx, dy, gx, gy);
    float a = aspect_rad(gx, gy) * (180.0f / 3.14159265358979323846f);
    out[gid] = fmod(a + 360.0f, 360.0f);
}
"""

_NORMALS_SRC = _GRAD_MSL + r"""
kernel void normals_map(device const float* dem    [[buffer(0)]],
                        device float*       out    [[buffer(1)]],
                        device const float* params [[buffer(2)]],
                        uint gid [[thread_position_in_grid]]) {
    int H = (int)params[0], W = (int)params[1];
    float dx = params[2], dy = params[3];
    uint npix = (uint)(H*W);
    if (gid >= npix) return;
    int iy = (int)(gid / (uint)W), ix = (int)(gid % (uint)W);
    float gx, gy; grad(dem, H, W, iy, ix, dx, dy, gx, gy);
    float slope_r = atan(sqrt(gx*gx + gy*gy));
    float aspect_r = aspect_rad(gx, gy);
    float e = sin(slope_r) * sin(aspect_r);
    float n = sin(slope_r) * cos(aspect_r);
    float u = cos(slope_r);
    float norm = sqrt(e*e + n*n + u*u);
    if (norm > 0.0f) { e /= norm; n /= norm; u /= norm; }
    else             { e = 0.0f; n = 0.0f; u = 0.0f; }
    out[0*npix + gid] = e;
    out[1*npix + gid] = n;
    out[2*npix + gid] = u;
}
"""

_KERNELS = {}


def _raw_kernel(src, name):
    kern = _KERNELS.get(name)
    if kern is None:
        import macmetalpy as xp
        kern = xp.RawKernel(src, name)
        _KERNELS[name] = kern
    return kern


def _metal_grad_op(dem, dx, dy, src, name, bands=1):
    """Launch a per-pixel (slope/aspect/normals) kernel; one thread per pixel.
    ``bands`` is the number of output planes (3 for normals, else 1)."""
    import macmetalpy as xp

    H, W = dem.shape
    kern = _raw_kernel(src, name)
    d_dem = xp.asarray(np.ascontiguousarray(dem, np.float32))
    shape = (bands, H, W) if bands > 1 else (H, W)
    d_out = xp.empty(shape, dtype=np.float32)
    d_par = xp.asarray(np.array([H, W, dx, dy], np.float32))
    kern(H * W, [d_dem, d_out, d_par])
    return d_out.get().reshape(shape)


def metal_slope(dem, dx, dy):
    return _metal_grad_op(dem, dx, dy, _SLOPE_SRC, "slope_map")


def metal_aspect(dem, dx, dy):
    return _metal_grad_op(dem, dx, dy, _ASPECT_SRC, "aspect_map")


def metal_normals(dem, dx, dy):
    return _metal_grad_op(dem, dx, dy, _NORMALS_SRC, "normals_map", bands=3)
