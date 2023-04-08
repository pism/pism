#ifndef PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H
#define PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H

#include "pism/util/connected_components_lakecc.hh"

namespace pism {

class IceModelVec2CellType;

/*!
 * This class determines the maximum fill height of the lake basins by iteratively
 * checking the entire domain for a set of increasing water levels, as described in Hinck
 * et al. (2020)
 */
class LakeLevelCC : public FillingAlgCC<ValidCC<SinkCC> > {
public:
  LakeLevelCC(IceGrid::ConstPtr g, double density_ratio, const IceModelVec2S &bed,
              const IceModelVec2S &thk, const IceModelVec2Int &pism_mask);

  LakeLevelCC(IceGrid::ConstPtr g, double density_ratio, const IceModelVec2S &bed,
              const IceModelVec2S &thk, const IceModelVec2Int &pism_mask,
              const IceModelVec2Int &valid_mask);

  virtual ~LakeLevelCC() = default;

  inline void computeLakeLevel(double zMin, double zMax, double dz, IceModelVec2S &result) {
    computeLakeLevel(zMin, zMax, dz, connected_components::invalid, result);
  }

private:
  void prepare_mask(const IceModelVec2CellType &pism_mask, IceModelVec2Int &result);

  void labelMap(int run_number, const VecList &lists, IceModelVec2S &result) const;

  void fill2Level(double level, IceModelVec2S &result);

  bool ForegroundCond(int i, int j) const;

  void computeLakeLevel(double zMin, double zMax, double dz, double offset, IceModelVec2S &result);

  double m_offset, m_level;
};

} // namespace pism

#endif /* PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H */
