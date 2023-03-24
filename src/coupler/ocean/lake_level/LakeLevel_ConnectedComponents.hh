#ifndef _LAKELEVEL_CONNECTEDCOMPONENTS_H_
#define _LAKELEVEL_CONNECTEDCOMPONENTS_H_

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


/*!
 * This class marks cells as invalid that are either ice-covered or not connected by an
 * ice-free corridor within the domain margin. It is used to restrict the formation of
 * lakes in the ice sheet interior and of subglacial lakes where thin ice covers a deep
 * basin.
 */
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


/*!
 * To remove narrow lakes, which are often related to under-resolved topography, the
 * FilterLakesCC class checks the lakes’ geometry. Only lakes containing one (or more)
 * cells that have at least a certain amount of neighbors that also are part of that lake
 * are retained
 */
class FilterLakesCC : public ValidCC<ConnectedComponents> {
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


/*!
 * LakePropertiesCC collects the minimum and maximum current water level of each lake
 * basin.
 */
class LakePropertiesCC : public ConnectedComponents {
public:
  LakePropertiesCC(IceGrid::ConstPtr g, const double fill_value, const IceModelVec2S &target_level,
                   const IceModelVec2S &lake_level);
  ~LakePropertiesCC();
  void getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level);

private:
  const double m_fill_value;
  const IceModelVec2S *m_target_level, *m_current_level;
  IceModelVec2S m_min_lakelevel, m_max_lakelevel;

  void setRunMinLevel(double level, int run, VecList &lists);
  void setRunMaxLevel(double level, int run, VecList &lists);
  inline bool isLake(const double level) const {
    return (level != m_fill_value);
  }

  inline bool ForegroundCond(const double target, const double current) const {
    return (isLake(target));
  }

protected:
  virtual void init_VecList(VecList &lists, const unsigned int length);
  virtual bool ForegroundCond(const int i, const int j) const;
  virtual void labelMask(int run_number, const VecList &lists);
  virtual void treatInnerMargin(const int i, const int j,
                                const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                VecList &lists, bool &changed);
  virtual void startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists);
  virtual void continueRun(const int i, const int j, int &run_number, VecList &lists);
};

class LakeAccumulatorCCSerial : public ConnectedComponentsSerial {
public:
  LakeAccumulatorCCSerial(IceGrid::ConstPtr g, const double fill_value);
  ~LakeAccumulatorCCSerial();
  void init(const IceModelVec2S &lake_level);
  void accumulate(const IceModelVec2S &in, IceModelVec2S &result);

protected:
  virtual bool ForegroundCond(const int i, const int j) const;

private:
  bool m_initialized;
  double m_fill_value;
  VecList m_lists;
  int m_run_number;

  void prepare_mask(const IceModelVec2S &lake_level);

  inline bool ForegroundCond(const int mask) const {
    return (mask > 0);
  }

  inline bool isLake(const double level) {
    return (level != m_fill_value);
  }
};

} // namespace pism

#endif
