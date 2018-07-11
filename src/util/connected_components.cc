#include "pism/util/pism_utilities.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/error_handling.hh"

#include "connected_components.hh"

namespace pism {

ConnectedComponentsBase::ConnectedComponentsBase(const int dList):
  m_dList(dList) {
  //empty
}

ConnectedComponentsBase::~ConnectedComponentsBase() {
  //empty
}

inline void ConnectedComponentsBase::resizeLists(VecList &lists, const int new_length) {
  for (VecList::iterator it = lists.begin(); it != lists.end(); it++) {
    it->second.resize(new_length);
  }
}

inline void ConnectedComponentsBase::run_union(RunVec &parents, int run1, int run2) {
  if ((parents[run1] == run2) or (parents[run2] == run1)) {
    return;
  }

  run1 = trackParentRun(run1, parents);
  run2 = trackParentRun(run2, parents);

  if (run1 > run2) {
    parents[run1] = run2;
  } else if (run1 < run2) {
    parents[run2] = run1;
  }
}

inline void ConnectedComponentsBase::check_cell(const int i, const int j, const bool isWest, const bool isSouth, const int mask_w, const int mask_s, int &run_number, VecList &lists, unsigned int &max_items) {
  //Check Foreground Pixel
  if (not isWest and (mask_w > 0)) {
    // west neighbor is also foreground: continue the run
    continueRun(i, j, run_number, lists);
  } else {
    //west neighbor is a background pixel (or this is westmost column): start a new run
    int parent;
    if (not isSouth and (mask_s > 0)) {
      //check the pixel south and set the parent
      parent = mask_s;
    } else {
      parent = 0;
    }
    startNewRun(i, j, run_number, parent, lists);
  }

  if (not isSouth and (mask_s > 0)) {
    mergeRuns(run_number, mask_s, lists);
  }

  //resize vectors if 'max_items' are exceeded
  if ((run_number + 1) >= max_items) {
    max_items += m_dList;
    resizeLists(lists, max_items);
  }
}

int ConnectedComponentsBase::trackParentRun(int run, const RunVec &parents) {
  while (parents[run] != 0) {
    run = parents[run];
  }
  return run;
}

void ConnectedComponentsBase::init_VecList(VecList &lists, const unsigned int size) {
  RunVec parents(size), lengths(size), j_vec(size), i_vec(size);
  lists["parents"] = parents;
  lists["lengths"] = lengths;
  lists["j"] = j_vec;
  lists["i"] = i_vec;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["parents"][k] = 0;
    lists["lengths"][k] = 0;
    lists["j"][k] = 0;
    lists["i"][k] = 0;
  }
}

void ConnectedComponentsBase::startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists) {
  run_number += 1;
  lists["i"][run_number] = i;
  lists["j"][run_number] = j;
  lists["lengths"][run_number] = 1;
  lists["parents"][run_number] = parent;
}

void ConnectedComponentsBase::continueRun(const int i, const int j, int &run_number, VecList &lists) {
  lists["lengths"][run_number] += 1;
}

void ConnectedComponentsBase::mergeRuns(const int run_number, const int run_south, VecList &lists) {
  run_union(lists["parents"], run_south, run_number);
}



ConnectedComponents::ConnectedComponents(IceGrid::ConstPtr g)
    : m_grid(g),
      ConnectedComponentsBase(m_grid->ym()),
      m_i_local_first(m_grid->xs()),
      m_i_local_last(m_i_local_first + m_grid->xm() - 1),
      m_j_local_first(m_grid->ys()),
      m_j_local_last(m_j_local_first + m_grid->ym() - 1),
      m_i_global_first(0),
      m_i_global_last(m_grid->Mx() - 1),
      m_j_global_first(0),
      m_j_global_last(m_grid->My() - 1) {

  m_mask_run.create(m_grid, "mask_run", WITH_GHOSTS, 1);
  m_masks.push_back(&m_mask_run);
}

ConnectedComponents::~ConnectedComponents() {
  //empty
}

void ConnectedComponents::compute_runs(int &run_number, VecList &lists, unsigned int &max_items) {
  IceModelVec::AccessList list;
  addFieldVecAccessList(m_masks, list);
  addFieldVecAccessList(m_fields, list);

  //Assign Pixels to runs
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ForegroundCond(i, j)) {
      bool isWest = (i <= m_i_local_first), isSouth = (j <= m_j_local_first);
      const int mask_w = isWest  ? 0 : m_mask_run(i-1, j),
                mask_s = isSouth ? 0 : m_mask_run(i, j-1);
      check_cell(i, j, isWest, isSouth, mask_w, mask_s, run_number, lists, max_items);
      m_mask_run(i, j) = run_number;
    }
  }

  labelMask(run_number, lists);

  //iteratively adapt fields amongst processor domains
  bool updateAtBoundaries = true;
  while (updateAtBoundaries) {
    updateAtBoundaries = updateRunsAtBoundaries(lists);
    labelMask(run_number, lists);
  }
}

void ConnectedComponents::labelMask(int run_number, const VecList &lists) {
  IceModelVec::AccessList list;
  addFieldVecAccessList(m_masks, list);

  const RunVec &i_vec   = lists.find("i")->second,
               &j_vec   = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    for (unsigned int n = 0; n < len_vec[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      m_mask_run(i, j) = label;
    }
  }
}

bool ConnectedComponents::updateRunsAtBoundaries(VecList &lists) {
  IceModelVec::AccessList list;
  addFieldVecAccessList(m_masks, list);
  updateGhosts(m_masks);

  bool changed = false;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest   = ((i == m_i_local_first) and not(i == m_i_global_first)),
         isEast   = ((i == m_i_local_last) and not(i == m_i_global_last)),
         isSouth  = ((j == m_j_local_first) and not(j == m_j_global_first)),
         isNorth  = ((j == m_j_local_last) and not(j == m_j_global_last)),
         isMargin = (isWest or isEast or isSouth or isNorth);

    if (isMargin) {
      treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);
    }
  }

  return (GlobalOr(m_grid->com, changed));
}

void ConnectedComponents::addFieldVecAccessList(FieldVec &fields, IceModelVec::AccessList &list) {
  for (FieldVec::iterator it = fields.begin(); it != fields.end(); it++) {
    IceModelVec *field = *it;
    list.add(*field);
  }
}

void ConnectedComponents::addFieldVecAccessList(ConstFieldVec &fields, IceModelVec::AccessList &list) {
  for (ConstFieldVec::iterator it = fields.begin(); it != fields.end(); it++) {
    const IceModelVec *field = *it;
    list.add(*field);
  }
}

void ConnectedComponents::updateGhosts(FieldVec &in) {
  for (FieldVec::iterator it = in.begin(); it != in.end(); it++) {
    IceModelVec *field = *it;
    field->update_ghosts();
  }
}


ConnectedComponentsSerial::ConnectedComponentsSerial(int Mx, int My)
  :ConnectedComponentsBase(My),
   m_Mx(Mx),
   m_My(My) {

  // memory allocation
  PetscErrorCode ierr = 0;

  // m_mask_run
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_mask_run_vec.rawptr());;
  PISM_CHK(ierr, "VecCreateSeq");

  m_mask_run.reset(new petsc::VecArray2D(m_mask_run_vec, m_Mx, m_My));
}

ConnectedComponentsSerial::~ConnectedComponentsSerial() {
  //empty
}


void ConnectedComponentsSerial::compute_runs(int &run_number, VecList &lists, unsigned int &max_items) {

  for (int j = 0; j < m_My; j++) {
    for (int i = 0; i < m_Mx; i++) {
      if (ForegroundCond(i, j)) {
        bool isWest = (i <= 0), isSouth = (j <= 0);
        const int mask_w = isWest  ? 0 : (*m_mask_run)(i-1, j),
                  mask_s = isSouth ? 0 : (*m_mask_run)(i, j-1);
        check_cell(i, j, isWest, isSouth, mask_w, mask_s, run_number, lists, max_items);
        (*m_mask_run)(i, j) = run_number;
      }
    }
  }
  //I think this is not needed in the serial case. Might be run for further application
  //labelMask(run_number, lists);
}

void ConnectedComponentsSerial::labelMask(int run_number, const VecList &lists) {
  const RunVec &i_vec   = lists.find("i")->second,
               &j_vec   = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    for (unsigned int n = 0; n < len_vec[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      (*m_mask_run)(i, j) = label;
    }
  }
}


SinkCC::SinkCC(IceGrid::ConstPtr g)
  :ConnectedComponents(g) {
  //empty
}

SinkCC::~SinkCC() {
  //empty
}

void SinkCC::setRunSink(int run, RunVec &parents) {
  if ((run == 0) or (run == 1)) {
    return;
  }

  run = trackParentRun(run, parents);
  if (run != 1) {
    parents[run] = 1;
  }
}

bool SinkCC::SinkCond(const int i, const int j) {
  const int mask = m_mask_run(i, j);
  return (mask == 1);
}

void SinkCC::treatInnerMargin(const int i, const int j,
                              const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                              VecList &lists, bool &changed) {
  ConnectedComponents::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  const int run = m_mask_run.as_int(i, j);
  if (run > 1) {
    //Lake at inner boundary
    StarStencil<int> mask_star = m_mask_run.int_star(i, j);
    bool WestSink = (isWest and (mask_star.w == 1)), EastSink = (isEast and (mask_star.e == 1)),
         SouthSink = (isSouth and (mask_star.s == 1)), NorthSink = (isNorth and (mask_star.n == 1));

    if (WestSink or EastSink or SouthSink or NorthSink) {
      //Lake on other side overflowing
      lists["parents"][run] = 1;
      changed = true;
    }
  }
}

void SinkCC::startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists) {
  if (SinkCond(i, j)) {
    parent = 1;
  }
  ConnectedComponents::startNewRun(i, j, run_number, parent, lists);
}

void SinkCC::continueRun(const int i, const int j, int &run_number, VecList &lists) {
  ConnectedComponents::continueRun(i, j, run_number, lists);
  if (SinkCond(i, j)) {
    setRunSink(run_number, lists["parents"]);
  }
}



MaskCC::MaskCC(IceGrid::ConstPtr g)
  :SinkCC(g) {
  //empty
}

MaskCC::~MaskCC() {
  //empty
}

void MaskCC::compute_mask(IceModelVec2Int &mask) {
  m_mask_run.copy_from(mask);

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelMask(run_number, lists, mask);
}

bool MaskCC::ForegroundCond(const int i, const int j) const {
  const int mask = m_mask_run.as_int(i, j);
  return (mask > 0);
}

void MaskCC::labelMask(const int run_number, const VecList &lists, IceModelVec2Int &result) {
  IceModelVec::AccessList list{&result};
  result.set(0);

  const RunVec &i_vec = lists.find("i")->second,
               &j_vec = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    if (label > 1) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = 1;
      }
    }
  }
}

} //namespace pism