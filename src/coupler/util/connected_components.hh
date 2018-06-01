#ifndef _CONNECTEDCOMPONENTS_H_
#define _CONNECTEDCOMPONENTS_H_

#include <vector>
#include <map>

#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class ConnectedComponents {
public:
  ConnectedComponents(IceGrid::ConstPtr g);
  ~ConnectedComponents();

  typedef std::vector<double> RunVec;
  typedef std::map<std::string, RunVec > VecList;

private:
  const int m_i_local_first, m_i_local_last, m_j_local_first, m_j_local_last, m_i_global_first, m_i_global_last,
            m_j_global_first, m_j_global_last;

  void checkForegroundPixel(const int i, const int j, int &run_number, VecList &lists);
  void resizeLists(VecList &lists, const int new_length);
  void run_union(RunVec &parents, int run1, int run2);

protected:
  const IceGrid::ConstPtr m_grid;
  IceModelVec2Int m_mask_run;
  std::vector<const PetscAccessible*> m_masks, m_fields;

  void compute_runs(int &run_number, VecList &lists, unsigned int &max_items);
  virtual void init_VecList(VecList &lists, const unsigned int length);
  virtual bool ForegroundCond(const int i, const int j) const = 0;
  virtual void treatInnerMargin(const int i, const int j,
                                const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                VecList &lists, bool &changed) = 0;
  int trackParentRun(int run, const RunVec &parents);
  virtual void startNewRun(const int i, const int j, int &run_number, VecList &lists, const int parent);
  virtual void continueRun(const int i, const int j, int &run_number, VecList &lists);
  virtual void mergeRuns(const int run_number, const int run_south, VecList &lists);
  virtual void labelMask(int run_number, const VecList &lists);
  virtual bool updateRunsAtBoundaries(VecList &lists);
};

} //namespace pism

#endif
