// Standalone driver to exercise the real PISM terrain horizon kernel against the
// solshade-produced reference file. Compile with the kernel source (no PISM/PETSc needed):
//
//   c++ -O2 -std=c++11 -I <build>/include verify_horizon_main.cc \
//       ../../src/coupler/surface/terrain_insolation_kernel.cc -o verify_horizon
//
// Reads a raw float64 DEM (row-major dem[j*Mx + i], j increasing north) and writes the
// horizon (radians) for a sub-rectangle of cells as raw float64, shape
// (n_dir, j1-j0, i1-i0). Azimuth direction k is 2*pi*k/n_dir, clockwise from north.

#include "pism/coupler/surface/terrain_insolation_kernel.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace pism::surface::terrain;

int main(int argc, char **argv) {
  if (argc != 13) {
    std::fprintf(stderr,
                 "usage: %s dem.bin Mx My dx dy n_dir max_dist step i0 i1 j0 j1 > out.bin\n",
                 argv[0]);
    return 1;
  }
  const char *dem_path = argv[1];
  int Mx = std::atoi(argv[2]);
  int My = std::atoi(argv[3]);
  double dx = std::atof(argv[4]);
  double dy = std::atof(argv[5]);
  int n_dir = std::atoi(argv[6]);
  double max_dist = std::atof(argv[7]);
  double step = std::atof(argv[8]);
  int i0 = std::atoi(argv[9]), i1 = std::atoi(argv[10]);
  int j0 = std::atoi(argv[11]), j1 = std::atoi(argv[12]);

  std::vector<double> dem((size_t)Mx * My);
  FILE *f = std::fopen(dem_path, "rb");
  if (!f || std::fread(dem.data(), sizeof(double), dem.size(), f) != dem.size()) {
    std::fprintf(stderr, "failed to read DEM\n");
    return 1;
  }
  std::fclose(f);

  std::vector<double> az(n_dir);
  for (int k = 0; k < n_dir; ++k) {
    az[k] = 2.0 * M_PI * k / n_dir;
  }

  // output order (k, j, i)
  for (int k = 0; k < n_dir; ++k) {
    for (int j = j0; j < j1; ++j) {
      for (int i = i0; i < i1; ++i) {
        double h = ray_horizon(dem.data(), Mx, My, dx, dy, i, j, az[k], step, max_dist);
        std::fwrite(&h, sizeof(double), 1, stdout);
      }
    }
  }
  return 0;
}
