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
  list.add(m_masks);
  list.add(m_fields);

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

  bool updateAtBoundaries = true;

  //iteratively adapt fields amongst processor domains
  while (updateAtBoundaries) {
    updateAtBoundaries = updateRunsAtBoundaries(lists);

    labelMask(run_number, lists);
  }
}

void ConnectedComponents::init_VecList(VecList &lists, const unsigned int size) {
  RunVec parents(size), lengths(size), j_vec(size), i_vec(size);
  lists["parents"] = parents;
  lists["lengths"] = lengths;
  lists["j_vec"]   = j_vec;
  lists["i_vec"]   = i_vec;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["parents"][k] = 0;
    lists["lengths"][k] = 0;
    lists["j_vec"][k]   = 0;
    lists["i_vec"][k]   = 0;
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

    startNewRun(i, j, run_number, lists, parent);
  }

  if (not isSouth and (mask_star.s > 0)) {
    mergeRuns(run_number, mask_star.s, lists);
  }

}

void ConnectedComponents::startNewRun(const int i, const int j, int &run_number, VecList &lists, const int parent) {
  run_number += 1;
  lists["i_vec"][run_number]   = i;
  lists["j_vec"][run_number]   = j;
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
  list.add(m_masks);

  const RunVec i_vec = lists.find("i_vec")->second,
               j_vec = lists.find("j_vec")->second;

  for (int k = 0; k <= run_number; ++k) {
    int label = trackParentRun(k, lists.find("parents")->second);
    for (unsigned int n = 0; n < lists.find("lengths")->second[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      m_mask_run(i, j) = label;
    }
  }
}

bool ConnectedComponents::updateRunsAtBoundaries(VecList &lists) {
  (void) lists;

  return false;
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


} //namespace pism