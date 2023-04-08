#ifndef _SEALEVEL_CONNECTEDCOMPONENTS_H_
#define _SEALEVEL_CONNECTEDCOMPONENTS_H_

#include "pism/util/connected_components_lakecc.hh"

namespace pism {

/*!
 * This class identifies basins in the topography that are below the global sea level but
 * are not connected to either of the domain margins.
 */
class SeaLevelCC : public FillingAlgCC<SinkCC> {
public:
  SeaLevelCC(IceGrid::ConstPtr g,
             double density_ratio,
             const IceModelVec2S &bed,
             const IceModelVec2S &thk);
  virtual ~SeaLevelCC() = default;

  void computeMask(const IceModelVec2S &SeaLevel, double Offset, IceModelVec2Int &result);

protected:
  const IceModelVec2S *m_sea_level;
  double m_offset;

  virtual bool ForegroundCond(int i, int j) const;
};

} // end of the namespace pism

#endif
