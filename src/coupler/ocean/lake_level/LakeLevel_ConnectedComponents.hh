#ifndef _LAKELEVEL_CONNECTEDCOMPONENTS_H_
#define _LAKELEVEL_CONNECTEDCOMPONENTS_H_

#include "pism/util/IceModelVec2CellType.hh"
#include "pism/coupler/util/connected_components.hh"

namespace pism {

class LakeLevelCC : public FillingAlgCC<ValidSinkCC> {
public:
  LakeLevelCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed,
              const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, const double fill_value);
  LakeLevelCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed,
              const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, const double fill_value,
              const IceModelVec2Int &valid_mask);
  ~LakeLevelCC();
  void computeLakeLevel(const double zMin, const double zMax, const double dz, const double offset, IceModelVec2S &result);
  inline void computeLakeLevel(const double zMin, const double zMax, const double dz, IceModelVec2S &result) {
    computeLakeLevel(zMin, zMax, dz, m_fill_value, result);
  }

protected:
  void prepare_mask(const IceModelVec2CellType &pism_mask);
  void labelMap(const int run_number, const VecList &lists, IceModelVec2S &result);
  void fill2Level(const double level, IceModelVec2S &result);
  virtual bool ForegroundCond(const int i, const int j) const;

private:
  double m_offset, m_level;
};


class IsolationCC : public SinkCC {
public:
  IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk,
              const double thk_theshold);
  ~IsolationCC();
  void find_isolated_spots(IceModelVec2Int &result);

protected:
  virtual bool ForegroundCond(const int i, const int j) const;
  void labelIsolatedSpots(const int run_number, const VecList &lists, IceModelVec2Int &result);
  void prepare_mask();

private:
  const double m_thk_threshold;
  const IceModelVec2S *m_thk;
  inline bool ForegroundCond(const double thk, const int mask) const {
    return ((thk < m_thk_threshold) or (mask > 0));
  }
};


class FilterLakesCC : public ValidSinkCC {
public:
  FilterLakesCC(IceGrid::ConstPtr g, const double fill_value);
  ~FilterLakesCC();
  void filter_map(const int n_filter, IceModelVec2S &lake_level);

protected:
  virtual bool ForegroundCond(const int i, const int j) const;

private:
  double m_fill_value;

  void labelMap(const int run_number, const VecList &lists, IceModelVec2S &result);
  void prepare_mask(const IceModelVec2S &lake_level);
  void set_mask_validity(const int n_filter, const IceModelVec2S &lake_level);

  inline bool ForegroundCond(const int mask) const {
    return (mask > 1);
  }

  inline bool isLake(const double level) {
    return (level != m_fill_value);
  }
};

} // namespace pism

#endif
