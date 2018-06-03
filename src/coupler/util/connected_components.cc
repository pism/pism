#include "pism/util/pism_utilities.hh"
#include "pism/util/iceModelVec.hh"

#include "connected_components.hh"

namespace pism {

ConnectedComponents::ConnectedComponents(IceGrid::ConstPtr g)
    : m_grid(g),
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

      checkForegroundPixel(i, j, run_number, lists);

      m_mask_run(i, j) = run_number;

      //resize vectors if 'max_items' are exceeded
      if ((run_number + 1) >= max_items) {
        max_items += m_grid->ym();
        resizeLists(lists, max_items);
      }
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

void ConnectedComponents::init_VecList(VecList &lists, const unsigned int size) {
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

void ConnectedComponents::checkForegroundPixel(const int i, const int j, int &run_number, VecList &lists) {
  bool isWest = (i <= m_i_local_first), isSouth = (j <= m_j_local_first);
  StarStencil<int> mask_star = m_mask_run.int_star(i, j);

  if (not isWest and (mask_star.w > 0)) {
    // west neighbor is also foreground: continue the run
    continueRun(i, j, run_number, lists);
  } else {
    //west neighbor is a background pixel (or this is westmost column): start a new run
    int parent;
    if (not isSouth and (mask_star.s > 0)) {
      //check the pixel south and set the parent
      parent = (int)mask_star.s;
    } else {
      parent = 0;
    }
    startNewRun(i, j, run_number, parent, lists);
  }

  if (not isSouth and (mask_star.s > 0)) {
    mergeRuns(run_number, mask_star.s, lists);
  }
}

void ConnectedComponents::startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists) {
  run_number += 1;
  lists["i"][run_number] = i;
  lists["j"][run_number] = j;
  lists["lengths"][run_number] = 1;
  lists["parents"][run_number] = parent;
}

void ConnectedComponents::continueRun(const int i, const int j, int &run_number, VecList &lists) {
  lists["lengths"][run_number] += 1;
}

void ConnectedComponents::mergeRuns(const int run_number, const int run_south, VecList &lists) {
  run_union(lists["parents"], run_south, run_number);
}

void ConnectedComponents::resizeLists(VecList &lists, const int new_length) {
  for (VecList::iterator it = lists.begin(); it != lists.end(); it++) {
    it->second.resize(new_length);
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
      treatInnerMargin(i, j, isNorth, isEast, isWest, isWest, lists, changed);
    }
  }

  return (GlobalOr(m_grid->com, changed));
}

void ConnectedComponents::run_union(RunVec &parents, int run1, int run2) {
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

int ConnectedComponents::trackParentRun(int run, const RunVec &parents) {
  while (parents[run] != 0) {
    run = parents[run];
  }
  return run;
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



ValidSinkCC::ValidSinkCC(IceGrid::ConstPtr g)
  :SinkCC(g) {
  m_mask_validity.create(m_grid, "mask_validity", WITH_GHOSTS, 1);
  m_masks.push_back(&m_mask_validity);
}

ValidSinkCC::~ValidSinkCC() {
  //empty
}

void ValidSinkCC::init_VecList(VecList &lists, const unsigned int size) {
  SinkCC::init_VecList(lists, size);

  RunVec valid_list(size);
  lists["valid"] = valid_list;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["valid"][k] = 1;
  }
}

void ValidSinkCC::setRunValid(int run, VecList &lists) {
  if ((run == 0) or (run == 1)) {
    return;
  }

  run = trackParentRun(run, lists["parents"]);
  if (run != 1) {
    lists["valid"][run] = 1;
  }
}

void ValidSinkCC::treatInnerMargin(const int i, const int j,
                                   const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                   VecList &lists, bool &changed) {
  SinkCC::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  const int run = m_mask_run.as_int(i, j);
  if (run > 1) {
    //Lake at inner boundary
    const bool isValid = (lists["valid"][run] > 0);
    if (not isValid) {
      //Lake at this side is not labeled as valid
      StarStencil<int> mask_isValid_star = m_mask_validity.int_star(i, j);

      bool WestValid  = (isWest  and (mask_isValid_star.w == 1)),
           EastValid  = (isEast  and (mask_isValid_star.e == 1)),
           SouthValid = (isSouth and (mask_isValid_star.s == 1)),
           NorthValid = (isNorth and (mask_isValid_star.n == 1));

      if (WestValid or EastValid or SouthValid or NorthValid) {
        //Lake at other side is not completely covered with ice
        lists["valid"][run] = 1;
        changed = true;
      }
    }
  }
}

void ValidSinkCC::startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists) {
  SinkCC::startNewRun(i, j, run_number, parent, lists);
  const bool isValid = (m_mask_validity(i, j) > 0);
  lists["valid"][run_number] = isValid ? 1 : 0;
}

void ValidSinkCC::continueRun(const int i, const int j, int &run_number, VecList &lists) {
  SinkCC::continueRun(i, j, run_number, lists);
  const bool isValid = (m_mask_validity(i, j) > 0);
  if (isValid) {
    setRunValid(run_number, lists);
  }
}

void ValidSinkCC::mergeRuns(const int run_number, const int run_south, VecList &lists) {
  SinkCC::mergeRuns(run_number, run_south, lists);
  const bool isValid = (lists["valid"][run_number] > 0);
  if (isValid) {
    setRunValid(run_number, lists);
  }
}

void ValidSinkCC::labelMask(int run_number, const VecList &lists) {
  IceModelVec::AccessList list;
  addFieldVecAccessList(m_masks, list);

  const RunVec &i_vec = lists.find("i")->second,
               &j_vec = lists.find("j")->second,
               &len_vec   = lists.find("lengths")->second,
               &parents   = lists.find("parents")->second,
               &valid_vec = lists.find("valid")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    const int label_valid = valid_vec[label];
    for (unsigned int n = 0; n < len_vec[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      m_mask_run(i, j) = label;
      m_mask_validity(i, j) = label_valid;
    }
  }
}

} //namespace pism