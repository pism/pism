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
  typedef std::vector<IceModelVec*> FieldVec;
  typedef std::vector<const IceModelVec*> ConstFieldVec;

private:
  void checkForegroundPixel(const int i, const int j, int &run_number, VecList &lists);
  void resizeLists(VecList &lists, const int new_length);
  void run_union(RunVec &parents, int run1, int run2);

protected:
  const IceGrid::ConstPtr m_grid;
  IceModelVec2Int m_mask_run;
  FieldVec m_masks;
  ConstFieldVec m_fields;
  const int m_i_local_first, m_i_local_last, m_j_local_first, m_j_local_last,
            m_i_global_first, m_i_global_last, m_j_global_first, m_j_global_last;

  void compute_runs(int &run_number, VecList &lists, unsigned int &max_items);
  bool updateRunsAtBoundaries(VecList &lists);
  int trackParentRun(int run, const RunVec &parents);
  void updateGhosts(FieldVec &in);
  void addFieldVecAccessList(FieldVec &field, IceModelVec::AccessList &list);
  void addFieldVecAccessList(ConstFieldVec &field, IceModelVec::AccessList &list);
  virtual void init_VecList(VecList &lists, const unsigned int length);
  virtual void startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists);
  virtual void continueRun(const int i, const int j, int &run_number, VecList &lists);
  virtual void mergeRuns(const int run_number, const int run_south, VecList &lists);
  virtual void labelMask(int run_number, const VecList &lists);
  virtual bool ForegroundCond(const int i, const int j) const = 0;
  virtual void treatInnerMargin(const int i, const int j,
                                const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                VecList &lists, bool &changed) = 0;
};



class SinkCC : public ConnectedComponents {
public:
  SinkCC(IceGrid::ConstPtr g);
  ~SinkCC();

private:
  void setRunSink(int run, RunVec &parents);

protected:
  virtual bool SinkCond(const int i, const int j);
  virtual void treatInnerMargin(const int i, const int j,
                                const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                VecList &lists, bool &changed);
  virtual void startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists);
  virtual void continueRun(const int i, const int j, int &run_number, VecList &lists);
};


class ValidSinkCC : public SinkCC {
public:
  ValidSinkCC(IceGrid::ConstPtr g);
  ~ValidSinkCC();

private:
  void setRunValid(int run, VecList &lists);

protected:
  IceModelVec2Int m_mask_validity;
  virtual void init_VecList(VecList &lists, const unsigned int length);
  virtual void treatInnerMargin(const int i, const int j,
                                const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                VecList &lists, bool &changed);
  virtual void startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists);
  virtual void continueRun(const int i, const int j, int &run_number, VecList &lists);
  virtual void mergeRuns(const int run_number, const int run_south, VecList &lists);
  virtual void labelMask(int run_number, const VecList &lists);
};



template<class CC>
class FillingAlgCC : public CC {
public:
  FillingAlgCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed, const IceModelVec2S &thk, const double fill_value);
  ~FillingAlgCC();
protected:
  double m_drho, m_fill_value;
  const IceModelVec2S *m_bed, *m_thk;
  virtual bool ForegroundCond(const double bed, const double thk, const int mask, const double Level, const double Offset) const;
};

template<class CC>
FillingAlgCC<CC>::FillingAlgCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed, const IceModelVec2S &thk, const double fill_value)
  : CC(g), m_drho(drho), m_bed(&bed), m_thk(&thk), m_fill_value(fill_value) {
  CC::m_fields.push_back(m_bed);
  CC::m_fields.push_back(m_thk);
}

template<class CC>
FillingAlgCC<CC>::~FillingAlgCC() {
  //empty
}

template<class CC>
bool FillingAlgCC<CC>::ForegroundCond(const double bed, const double thk, const int mask, const double Level, const double Offset) const {
  if (mask > 0) {
    return true;
  }

  if (Level == m_fill_value) {
    return true;
  }

  double level = Level;
  if (Offset != m_fill_value) {
    level += Offset;
  }

  return ((bed + (m_drho * thk)) < level);
}

} //namespace pism

#endif
