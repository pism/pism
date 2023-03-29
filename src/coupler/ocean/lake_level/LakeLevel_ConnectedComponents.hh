#ifndef PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H
#define PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H

#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/connected_components_lakecc.hh"

namespace pism {

/*!
 * This class determines the maximum fill height of the lake basins by iteratively
 * checking the entire domain for a set of increasing water levels, as described in Hinck
 * et al. (2020)
 */
class LakeLevelCC : public FillingAlgCC<ValidCC<SinkCC> > {
public:
  LakeLevelCC(IceGrid::ConstPtr g, double drho, const IceModelVec2S &bed,
              const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, double fill_value);
  LakeLevelCC(IceGrid::ConstPtr g, double drho, const IceModelVec2S &bed,
              const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, double fill_value,
              const IceModelVec2Int &valid_mask);
  ~LakeLevelCC();
  void computeLakeLevel(double zMin, double zMax, double dz, double offset, IceModelVec2S &result);
  inline void computeLakeLevel(double zMin, double zMax, double dz, IceModelVec2S &result) {
    computeLakeLevel(zMin, zMax, dz, m_fill_value, result);
  }

protected:
  void prepare_mask(const IceModelVec2CellType &pism_mask);
  void labelMap(int run_number, const VecList &lists, IceModelVec2S &result);
  void fill2Level(double level, IceModelVec2S &result);
  virtual bool ForegroundCond(int i, int j) const;

private:
  double m_offset, m_level;
};




class LakeAccumulatorCCSerial : public ConnectedComponentsSerial {
public:
  LakeAccumulatorCCSerial(IceGrid::ConstPtr g, double fill_value);
  virtual ~LakeAccumulatorCCSerial() = default;
  void init(const IceModelVec2S &lake_level);
  void accumulate(const IceModelVec2S &in, IceModelVec2S &result);

protected:
  virtual bool ForegroundCond(int i, int j) const;

private:
  bool m_initialized;
  double m_fill_value;
  VecList m_lists;
  int m_run_number;

  void prepare_mask(const IceModelVec2S &lake_level);

  inline bool isLake(double level) const {
    return (level != m_fill_value);
  }
};

} // namespace pism

#endif /* PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H */
