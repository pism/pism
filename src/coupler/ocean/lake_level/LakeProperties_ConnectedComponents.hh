#ifndef _LAKEPROPERTIESALG_CONNECTEDCOMPONENTS_H_
#define _LAKEPROPERTIESALG_CONNECTEDCOMPONENTS_H_

#include <vector>

#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class LakePropertiesCC {
public:
  LakePropertiesCC(IceGrid::ConstPtr g, const double drho, const double fill_value, const IceModelVec2S &target_level, const IceModelVec2S &lake_level, const IceModelVec2S &floating_thresh);
  ~LakePropertiesCC();
  void getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level, IceModelVec2S &min_float_level);

private:
  void checkForegroundPixel(unsigned int i, unsigned int j, StarStencil<int> mask_star, unsigned int &run_number, const double lake_level, const double float_thresh, std::vector<unsigned int> &i_vec, std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list);
  void labelMask(unsigned int run_number, std::vector<unsigned int> &rows, std::vector<unsigned int> &columns, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list);
//   void labelMap(unsigned int run_number, std::vector<unsigned int> &rows, std::vector<unsigned int> &columns, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths);
  bool updateRunsAtBoundaries(std::vector<unsigned int> &parents, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list);
  void run_union(std::vector<unsigned int> &parents, unsigned int run1, unsigned int run2);
  bool ForegroundCond(double target_ll, double current_ll) const;
  void setRunMinLevel(unsigned int run, std::vector<unsigned int> &parents, std::vector<double> &run_level_list, double level);
  void setRunMaxLevel(unsigned int run, std::vector<unsigned int> &parents, std::vector<double> &run_level_list, double level);

protected:
  double m_drho, m_fill_value;
  const IceGrid::ConstPtr m_grid;
  const int m_i_local_first, m_i_local_last, m_j_local_first, m_j_local_last, m_i_global_first, m_i_global_last,
      m_j_global_first, m_j_global_last;
  const IceModelVec2S *m_target_level, *m_current_level, *m_floating_threshold_level;
  IceModelVec2Int m_mask_run;
  IceModelVec2S m_min_ll_tmp, m_max_ll_tmp, m_min_float_tmp;

  unsigned int trackParentRun(unsigned int run, std::vector<unsigned int> &parents);

//   void labelMap_impl(unsigned int run_number, std::vector<unsigned int> &rows, std::vector<unsigned int> &columns, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths) = 0;
  void labelMask_impl(unsigned int run_number, std::vector<unsigned int> &rows, std::vector<unsigned int> &columns, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list);
  bool ForegroundCond_impl(double target_ll, double current_ll) const;
};

} //namespace pism

#endif
