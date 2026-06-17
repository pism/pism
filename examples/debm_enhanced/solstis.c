/* solstis.c -- terrain ops reproducing solshade, pure C / C-ABI, OpenMP.
 *
 *  - horizon_map : per-pixel horizon-angle map by self-contained ray-march
 *                  (solshade compute_horizon_map).
 *  - slope_map / aspect_map / normals_map : local-stencil slope, aspect and
 *                  ENU unit-normal vectors (solshade compute_slope_aspect_normals).
 *
 * Compiled to libsolstis.dylib by build.sh and called via ctypes.
 */
#include <math.h>
#include <stddef.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Bilinear sample of the DEM at fractional (row=fy, col=fx), replicating
 * scipy.ndimage.map_coordinates(order=1, mode='constant', cval=NaN):
 * any tap with nonzero weight that is out of bounds OR itself NaN makes the
 * whole sample NaN. */
static inline float bilinear(const float *dem, int H, int W, float fy, float fx)
{
    float ffy = floorf(fy), ffx = floorf(fx);
    int y0 = (int)ffy, x0 = (int)ffx;
    float ty = fy - ffy, tx = fx - ffx;

    float w00 = (1.0f - ty) * (1.0f - tx);
    float w01 = (1.0f - ty) * tx;
    float w10 = ty * (1.0f - tx);
    float w11 = ty * tx;

    float acc = 0.0f;

#define TAP(yy, xx, ww)                                                        \
    do {                                                                       \
        if ((ww) > 0.0f) {                                                     \
            if ((yy) < 0 || (yy) >= H || (xx) < 0 || (xx) >= W)                \
                return NAN;                                                     \
            float v = dem[(size_t)(yy) * (size_t)W + (size_t)(xx)];            \
            if (isnan(v))                                                      \
                return NAN;                                                     \
            acc += (ww) * v;                                                   \
        }                                                                      \
    } while (0)

    TAP(y0,     x0,     w00);
    TAP(y0,     x0 + 1, w01);
    TAP(y0 + 1, x0,     w10);
    TAP(y0 + 1, x0 + 1, w11);
#undef TAP

    return acc;
}

/* dem:        (H, W)        row-major float32
 * dx_pix,     (n_dir, ns)   per-direction, per-distance column offset (pixels)
 * dy_pix:     (n_dir, ns)   per-direction, per-distance row offset (pixels)
 * distances:  (ns,)         ray-sample distances in metres (distances[0] == 0)
 * out:        (n_dir, H, W) horizon angle in degrees, NaN where undefined
 */
void horizon_map(const float *dem, int H, int W,
                 const float *dx_pix, const float *dy_pix,
                 const float *distances, int n_dir, int ns,
                 float *out)
{
    const size_t plane = (size_t)H * (size_t)W;
    const float rad2deg = (float)(180.0 / M_PI);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {
            float elev0 = bilinear(dem, H, W, (float)iy, (float)ix);
            int bad0 = isnan(elev0);

            for (int idir = 0; idir < n_dir; ++idir) {
                float best = -INFINITY;
                if (!bad0) {
                    const float *dxr = dx_pix + (size_t)idir * ns;
                    const float *dyr = dy_pix + (size_t)idir * ns;
                    for (int k = 1; k < ns; ++k) {
                        float s = bilinear(dem, H, W,
                                           (float)iy + dyr[k],
                                           (float)ix + dxr[k]);
                        if (isnan(s))
                            continue;
                        float ang = atan2f(s - elev0, distances[k]);
                        if (ang > best)
                            best = ang;
                    }
                }
                out[(size_t)idir * plane + (size_t)iy * W + ix] =
                    isfinite(best) ? best * rad2deg : NAN;
            }
        }
    }
}


/* Elevation gradient at pixel (iy, ix), replicating numpy.gradient with the
 * default edge_order=1: central differences in the interior, first-order
 * one-sided differences on the first/last row/column. dx, dy are pixel sizes
 * (metres). dzdx is d/dx (columns), dzdy is d/dy (rows). */
static inline void grad(const float *dem, int H, int W, int iy, int ix,
                        float dx, float dy, float *dzdx, float *dzdy)
{
    const float *row = dem + (size_t)iy * W;
    float gx;
    if (W < 2)            gx = 0.0f;
    else if (ix == 0)     gx = (row[1] - row[0]) / dx;
    else if (ix == W - 1) gx = (row[W - 1] - row[W - 2]) / dx;
    else                  gx = (row[ix + 1] - row[ix - 1]) / (2.0f * dx);

    float gy;
    if (H < 2)            gy = 0.0f;
    else if (iy == 0)     gy = (dem[(size_t)W + ix]            - dem[ix]) / dy;
    else if (iy == H - 1) gy = (dem[(size_t)(H - 1) * W + ix]  - dem[(size_t)(H - 2) * W + ix]) / dy;
    else                  gy = (dem[(size_t)(iy + 1) * W + ix] - dem[(size_t)(iy - 1) * W + ix]) / (2.0f * dy);

    *dzdx = gx;
    *dzdy = gy;
}

/* Aspect angle in radians, atan2(-dzdx, dzdy), with the flat case (both
 * gradients exactly zero) pinned to 0 to match numpy's atan2(0, 0) == 0. */
static inline float aspect_rad(float gx, float gy)
{
    return (gx == 0.0f && gy == 0.0f) ? 0.0f : atan2f(-gx, gy);
}

/* Slope in degrees: angle from horizontal, atan(hypot(dzdx, dzdy)). */
void slope_map(const float *dem, int H, int W, float dx, float dy, float *out)
{
    const float rad2deg = (float)(180.0 / M_PI);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {
            float gx, gy;
            grad(dem, H, W, iy, ix, dx, dy, &gx, &gy);
            out[(size_t)iy * W + ix] = atanf(hypotf(gx, gy)) * rad2deg;
        }
    }
}

/* Aspect in degrees clockwise from North: (deg(atan2(-dzdx, dzdy)) + 360) % 360. */
void aspect_map(const float *dem, int H, int W, float dx, float dy, float *out)
{
    const float rad2deg = (float)(180.0 / M_PI);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {
            float gx, gy;
            grad(dem, H, W, iy, ix, dx, dy, &gx, &gy);
            float a = aspect_rad(gx, gy) * rad2deg;
            out[(size_t)iy * W + ix] = fmodf(a + 360.0f, 360.0f);
        }
    }
}

/* ENU unit normal vector, out shape (3, H, W): bands [east, north, up].
 *   E = sin(slope) sin(aspect), N = sin(slope) cos(aspect), U = cos(slope),
 * then L2-normalized (norm is ~1 analytically; matches solshade's defensive
 * normalize, with 0 where the norm is 0). */
void normals_map(const float *dem, int H, int W, float dx, float dy, float *out)
{
    const size_t plane = (size_t)H * (size_t)W;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int iy = 0; iy < H; ++iy) {
        for (int ix = 0; ix < W; ++ix) {
            float gx, gy;
            grad(dem, H, W, iy, ix, dx, dy, &gx, &gy);
            float slope_r = atanf(hypotf(gx, gy));
            float aspect_r = aspect_rad(gx, gy);
            float e = sinf(slope_r) * sinf(aspect_r);
            float n = sinf(slope_r) * cosf(aspect_r);
            float u = cosf(slope_r);
            float norm = sqrtf(e * e + n * n + u * u);
            if (norm > 0.0f) { e /= norm; n /= norm; u /= norm; }
            else             { e = n = u = 0.0f; }
            size_t gid = (size_t)iy * W + ix;
            out[0 * plane + gid] = e;
            out[1 * plane + gid] = n;
            out[2 * plane + gid] = u;
        }
    }
}
