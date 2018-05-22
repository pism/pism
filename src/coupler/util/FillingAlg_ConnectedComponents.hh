#ifndef _FILLINGALG_CONNECTEDCOMPONENTS_H_
#define _FILLINGALG_CONNECTEDCOMPONENTS_H_

#include <vector>

#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class FillingAlgCC {
public:
  FillingAlgCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed, const IceModelVec2S &thk,
               const bool checkValidity, const double fill_value);
  ~FillingAlgCC();
  void fill2Level();
  void fill2Level(double Level);
  void fill2Level(const IceModelVec2S &Level);
  void fill2Level(double Level, double Offset);
  void fill2Level(const IceModelVec2S &Level, double Offset);
  void fill2Level(const IceModelVec2S &Level, const IceModelVec2S &Offset);
  void setLevel(double Level);
  void setLevel(const IceModelVec2S &Level);
  void setOffset(double Offset);
  void setOffset(const IceModelVec2S &Offset);

private:
  void checkForegroundPixel(unsigned int i, unsigned int j, bool isSink, bool isValid, StarStencil<int> mask_star,
                            unsigned int &run_number, std::vector<unsigned int> &i_vec,
                            std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents,
                            std::vector<unsigned int> &lengths, std::vector<bool> &isValidList);
  void labelMask(unsigned int run_number, std::vector<unsigned int> &rows, std::vector<unsigned int> &columns,
                 std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<bool> &isValidList);
  void labelMap(unsigned int run_number, std::vector<unsigned int> &rows,
                std::vector<unsigned int> &columns, std::vector<unsigned int> &parents,
                std::vector<unsigned int> &lengths, std::vector<bool> &isValidList);
  bool updateRunsAtBoundaries(std::vector<unsigned int> &parents, std::vector<bool> &isValidList);
  void run_union(std::vector<unsigned int> &parents, unsigned int run1, unsigned int run2);
  void setRunSink(unsigned int run, std::vector<unsigned int> &parents);
  void setRunValid(unsigned int run, std::vector<unsigned int> &parents, std::vector<bool> &isValidList);
  bool SinkCond(int mask) const;
  bool ForegroundCond(double bed, double thk, int mask, double Level, double Offset) const;

  const bool m_check_validity;

protected:
  double m_drho, m_fill_value;
  const IceGrid::ConstPtr m_grid;
  const int m_i_local_first, m_i_local_last, m_j_local_first, m_j_local_last, m_i_global_first, m_i_global_last,
            m_j_global_first, m_j_global_last;
  IceModelVec2S m_floatation_level;
  const IceModelVec2S *m_bed, *m_thk;
  IceModelVec2Int m_mask_run, m_mask_validity, m_level, m_offset;

  void get_floatation_level(IceModelVec2S &result) const;
  void prepare_mask();
  unsigned int trackParentRun(unsigned int run, std::vector<unsigned int> &parents);

  virtual void prepare_mask_impl() = 0;
  virtual void labelMap_impl(unsigned int run_number, std::vector<unsigned int> &rows,
                             std::vector<unsigned int> &columns, std::vector<unsigned int> &parents,
                             std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) = 0;
  virtual void labelMask_impl(unsigned int run_number, std::vector<unsigned int> &rows,
                              std::vector<unsigned int> &columns, std::vector<unsigned int> &parents,
                              std::vector<unsigned int> &lengths, std::vector<bool> &isValidList);
  virtual bool SinkCond_impl(int mask) const;
  virtual bool ForegroundCond_impl(double bed, double thk, int mask, double Level, double Offset) const;
};

} //namespace pism

#endif
