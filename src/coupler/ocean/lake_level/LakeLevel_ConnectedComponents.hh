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


/*
class FilterLakesCC : public FillingAlgCC {
public:
  FilterLakesCC(IceGrid::ConstPtr g, const IceModelVec2S &lake_levels,
                const double fill_value);
  ~FilterLakesCC();
  void filter_map(int n_filter);
  void filtered_levels(IceModelVec2S &result) const;

protected:
  virtual bool ForegroundCond_impl(double lake_level, double thk, int mask,
                                   double Level, double Offset) const;
  virtual void labelMap_impl(unsigned int run_number,
                             std::vector<unsigned int> &i_vec,
                             std::vector<unsigned int> &j_vec,
                             std::vector<unsigned int> &parents,
                             std::vector<unsigned int> &lengths,
                             std::vector<bool> &isValidList);
  virtual void prepare_mask_impl();
  void set_mask_validity(int n_filter);
  const IceModelVec2S *m_lake_level;

};
*/
} // namespace pism

#endif
