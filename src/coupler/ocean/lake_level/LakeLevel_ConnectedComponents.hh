#ifndef _LAKELEVEL_CONNECTEDCOMPONENTS_H_
#define _LAKELEVEL_CONNECTEDCOMPONENTS_H_

#include <vector>

#include "pism/util/IceModelVec2CellType.hh"

#include "pism/coupler/util/FillingAlg_ConnectedComponents.hh"

namespace pism {

class LakeLevelCC : public FillingAlgCC {
public:
  LakeLevelCC(IceGrid::ConstPtr g, const double drho,
              const double thk_theshold, const IceModelVec2S &bed,
              const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
              IceModelVec2Int &run_mask, const double fill_value);
  LakeLevelCC(IceGrid::ConstPtr g, const double drho,
              const double thk_theshold, const IceModelVec2S &bed,
              const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
              const double fill_value);
  LakeLevelCC(IceGrid::ConstPtr g, const double drho,
              const double thk_theshold, const IceModelVec2S &bed,
              const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
              IceModelVec2Int &run_mask, const double fill_value,
              IceModelVec2Int &valid_mask);
  LakeLevelCC(IceGrid::ConstPtr g, const double drho,
              const double thk_theshold, const IceModelVec2S &bed,
              const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
              const double fill_value, IceModelVec2Int &valid_mask);
  ~LakeLevelCC();
  void floodMap(double zMin, double zMax, double dz);
  void lake_levels(IceModelVec2S &result) const;

protected:
  IceModelVec2CellType m_pism_mask;
  virtual void prepare_mask_impl();
  virtual void labelMap_impl(unsigned int run_number,
                             std::vector<unsigned int> &rows,
                             std::vector<unsigned int> &columns,
                             std::vector<unsigned int> &parents,
                             std::vector<unsigned int> &lengths,
                             std::vector<bool> &isValidList);
  void init(IceModelVec2Int &pism_mask);

private:
  const double m_thk_threshold;
};


class IsolationCC : public FillingAlgCC {
public:
  IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk,
              const double thk_theshold);
  ~IsolationCC();
  void find_isolated_spots();
  void isolation_mask(IceModelVec2Int &result) const;

protected:
  virtual bool ForegroundCond_impl(double bed, double thk, int mask,
                                   double Level, double Offset) const;
  virtual void labelMap_impl(unsigned int run_number,
                             std::vector<unsigned int> &i_vec,
                             std::vector<unsigned int> &j_vec,
                             std::vector<unsigned int> &parents,
                             std::vector<unsigned int> &lengths,
                             std::vector<bool> &isValidList);
  virtual void prepare_mask_impl();

private:
  const double m_thk_threshold;
};

} // namespace pism

#endif
