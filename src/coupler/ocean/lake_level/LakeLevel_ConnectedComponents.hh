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


/*!
 * This class marks cells as invalid that are either ice-covered or not connected by an
 * ice-free corridor within the domain margin. It is used to restrict the formation of
 * lakes in the ice sheet interior and of subglacial lakes where thin ice covers a deep
 * basin.
 */
class IsolationCC : public SinkCC {
public:
  IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk,
              double thk_theshold);
  ~IsolationCC();
  void find_isolated_spots(IceModelVec2Int &result);

protected:
  virtual bool ForegroundCond(int i, int j) const;
  void labelIsolatedSpots(int run_number, const VecList &lists, IceModelVec2Int &result);
  void prepare_mask();

private:
  const double m_thk_threshold;
  const IceModelVec2S *m_thk;

  inline bool ForegroundCond(double thk, int mask) const {
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
  FilterLakesCC(IceGrid::ConstPtr g, double fill_value);
  ~FilterLakesCC();
  void filter_map(int n_filter, IceModelVec2S &lake_level);

protected:
  virtual bool ForegroundCond(int i, int j) const;

private:
  double m_fill_value;

  void labelMap(int run_number, const VecList &lists, IceModelVec2S &result);
  void prepare_mask(const IceModelVec2S &lake_level);
  void set_mask_validity(int n_filter, const IceModelVec2S &lake_level);

  inline bool ForegroundCond(int mask) const {
    return (mask > 1);
  }

  inline bool isLake(double level) const {
    return (level != m_fill_value);
  }
};


/*!
 * LakePropertiesCC collects the minimum and maximum current water level of each lake
 * basin.
 */
class LakePropertiesCC : public ConnectedComponents {
public:
  LakePropertiesCC(IceGrid::ConstPtr g, double fill_value, const IceModelVec2S &target_level,
                   const IceModelVec2S &lake_level);
  ~LakePropertiesCC();
  void getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level);

private:
  const double m_fill_value;
  const IceModelVec2S *m_target_level, *m_current_level;
  IceModelVec2S m_min_lakelevel, m_max_lakelevel;

  void setRunMinLevel(double level, int run, VecList &lists);
  void setRunMaxLevel(double level, int run, VecList &lists);
  inline bool isLake(double level) const {
    return (level != m_fill_value);
  }

protected:
  virtual void init_VecList(VecList &lists, unsigned int length);
  virtual bool ForegroundCond(int i, int j) const;
  virtual void labelMask(int run_number, const VecList &lists);
  virtual void treatInnerMargin(int i, int j,
                                bool isNorth, bool isEast, bool isSouth, bool isWest,
                                VecList &lists, bool &changed);
  virtual void startNewRun(int i, int j, int &run_number, int &parent, VecList &lists);
  virtual void continueRun(int i, int j, int &run_number, VecList &lists);
};

class LakeAccumulatorCCSerial : public ConnectedComponentsSerial {
public:
  LakeAccumulatorCCSerial(IceGrid::ConstPtr g, double fill_value);
  ~LakeAccumulatorCCSerial();
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

  inline bool ForegroundCond(int mask) const {
    return (mask > 0);
  }

  inline bool isLake(double level) const {
    return (level != m_fill_value);
  }
};

} // namespace pism

#endif /* PISM_LAKE_LEVEL_CONNECTED_COMPONENTS_H */
