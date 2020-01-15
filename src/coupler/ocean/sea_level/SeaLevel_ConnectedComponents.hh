#ifndef _SEALEVEL_CONNECTEDCOMPONENTS_H_
#define _SEALEVEL_CONNECTEDCOMPONENTS_H_

#include "pism/util/connected_components_lakecc.hh"

namespace pism{

class SeaLevelCC : public FillingAlgCC<SinkCC> {
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
  void computeSeaLevel(IceModelVec2S &SeaLevel, const double Offset);
  void computeSeaLevel(const IceModelVec2S &SeaLevel, const double Offset, IceModelVec2S &result);
  void computeMask(const IceModelVec2S &SeaLevel, const double Offset, IceModelVec2Int &result);

protected:
  const IceModelVec2S *m_sea_level;
  double m_offset;
  void labelSLMap(const int run_number, const VecList lists, IceModelVec2S &result);
  void labelSLMask(const int run_number, const VecList lists, IceModelVec2Int &result);
  void prepare_mask();
  virtual bool ForegroundCond(const int i, const int j) const;
};

}

#endif
