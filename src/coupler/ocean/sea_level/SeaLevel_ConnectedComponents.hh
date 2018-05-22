#ifndef _SEALEVEL_CONNECTEDCOMPONENTS_H_
#define _SEALEVEL_CONNECTEDCOMPONENTS_H_

#include <vector>

#include "pism/coupler/util/FillingAlg_ConnectedComponents.hh"

namespace pism{

class SeaLevelCC : public FillingAlgCC {
public:
  SeaLevelCC(IceGrid::ConstPtr g,
             const double drho,
             const IceModelVec2S &bed,
             const IceModelVec2S &thk,
             const IceModelVec2Int &run_mask,
             const double fill_value);
  SeaLevelCC(IceGrid::ConstPtr g,
             const double drho,
             const IceModelVec2S &bed,
             const IceModelVec2S &thk,
             const double fill_value);
  ~SeaLevelCC();
  void fill2SeaLevel(double SeaLevel);
  void fill2SeaLevel(double SeaLevel, double Offset);
  void fill2SeaLevel(const IceModelVec2S &SeaLevel, double Offset);
  void fill2SeaLevel(const IceModelVec2S &SeaLevel, const IceModelVec2S &Offset);
  void sea_level_2D(IceModelVec2S &result) const;
  void sl_crop_mask(IceModelVec2Int& result) const;

protected:
  virtual void labelMap_impl(unsigned int run_number,
                             std::vector<unsigned int> &i_vec,
                             std::vector<unsigned int> &j_vec,
                             std::vector<unsigned int> &parents,
                             std::vector<unsigned int> &lengths,
                             std::vector<bool> &isValidList);
  virtual void prepare_mask_impl();
  IceModelVec2Int m_crop_mask;
};

}

#endif
