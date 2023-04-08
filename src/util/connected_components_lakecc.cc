
#include "Logger.hh"
#include "StarStencil.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/error_handling.hh"

#include "connected_components_lakecc.hh"
#include <functional>

namespace pism {

using connected_components::is_valid;

ConnectedComponentsBase::ConnectedComponentsBase(IceGrid::ConstPtr g):
  m_grid(g) {
  m_mask_run.create(m_grid, "mask_run", WITH_GHOSTS, 1);
}

inline void ConnectedComponentsBase::check_cell(int i, int j, bool isWest, bool isSouth, int mask_w, int mask_s, int &run_number, VecList &lists, unsigned int &max_items) {
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
  if ((run_number + 1) >= (int)max_items) {
    max_items += m_grid->ym();

    for (auto &pair : lists) {
      pair.second.resize(max_items);
    }
  }
}

void ConnectedComponentsBase::init_VecList(VecList &lists, unsigned int size) {
  std::vector<double> parents(size), lengths(size), j_vec(size), i_vec(size);
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

void ConnectedComponentsBase::startNewRun(int i, int j, int &run_number, int &parent, VecList &lists) {
  run_number += 1;
  lists["i"][run_number] = i;
  lists["j"][run_number] = j;
  lists["lengths"][run_number] = 1;
  lists["parents"][run_number] = parent;
}

void ConnectedComponentsBase::continueRun(int i, int j, int &run_number, VecList &lists) {
  (void) i;
  (void) j;
  lists["lengths"][run_number] += 1;
}

void ConnectedComponentsBase::mergeRuns(int run_number, int run_south, VecList &lists) {
  auto &parents = lists["parents"];
  int run1 = run_south;
  int run2 = run_number;
  {
    if ((parents[run1] == run2) or (parents[run2] == run1)) {
      return;
    }

    run1 = connected_components::trackParentRun(run1, parents);
    run2 = connected_components::trackParentRun(run2, parents);

    if (run1 > run2) {
      parents[run1] = run2;
    } else if (run1 < run2) {
      parents[run2] = run1;
    }
  }
}

ConnectedComponents::ConnectedComponents(IceGrid::ConstPtr g)
    : ConnectedComponentsBase(g),
      m_i_local_first(m_grid->xs()),
      m_i_local_last(m_i_local_first + m_grid->xm() - 1),
      m_j_local_first(m_grid->ys()),
      m_j_local_last(m_j_local_first + m_grid->ym() - 1),
      m_i_global_first(0),
      m_i_global_last(m_grid->Mx() - 1),
      m_j_global_first(0),
      m_j_global_last(m_grid->My() - 1) {

  m_masks.push_back(&m_mask_run);
}

void ConnectedComponents::compute_runs(int &run_number, VecList &lists, unsigned int &max_items) {
  IceModelVec::AccessList list;
  list.add(m_masks.begin(), m_masks.end());
  list.add(m_fields.begin(), m_fields.end());

  // Assign Pixels to runs
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ForegroundCond(i, j)) {
      bool
        isWest = (i <= m_i_local_first),
        isSouth = (j <= m_j_local_first);
      int
        mask_w = isWest  ? 0 : m_mask_run.as_int(i-1, j),
        mask_s = isSouth ? 0 : m_mask_run.as_int(i, j-1);
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
  connected_components::set_labels(run_number, lists, m_mask_run);
}

bool ConnectedComponents::updateRunsAtBoundaries(VecList &lists) {
  IceModelVec::AccessList list;
  list.add(m_masks.begin(), m_masks.end());

  for (auto *mask : m_masks) {
    mask->update_ghosts();
  }

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

ConnectedComponentsSerial::ConnectedComponentsSerial(IceGrid::ConstPtr g)
  : ConnectedComponentsBase(g) {
  m_mask_run_vec_p0 = m_mask_run.allocate_proc0_copy();
}

void ConnectedComponentsSerial::compute_runs(int &run_number, VecList &lists, unsigned int &max_items) {

  m_mask_run.put_on_proc0(*m_mask_run_vec_p0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      petsc::VecArray2D mask_run(*m_mask_run_vec_p0, m_grid->Mx(), m_grid->My());
      //We need a global pointer to the array to be able to access it from ForegroundCond(i, j)
      m_mask_run_p0_ptr = &mask_run;
      for (int j = 0; j < (int)m_grid->My(); j++) {
        for (int i = 0; i < (int)m_grid->Mx(); i++) {
          if (ForegroundCond(i, j)) {
            bool isWest = (i <= 0), isSouth = (j <= 0);
            const int mask_w = isWest  ? 0 : static_cast<int>(mask_run(i-1, j)),
                      mask_s = isSouth ? 0 : static_cast<int>(mask_run(i, j-1));
            check_cell(i, j, isWest, isSouth, mask_w, mask_s, run_number, lists, max_items);
            mask_run(i, j) = run_number;
          }
        }
      }
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}


SinkCC::SinkCC(IceGrid::ConstPtr g)
  :ConnectedComponents(g) {
  //empty
}

void SinkCC::setRunSink(int run, std::vector<double> &parents) {
  if ((run == 0) or (run == 1)) {
    return;
  }

  run = connected_components::trackParentRun(run, parents);
  if (run != 1) {
    parents[run] = 1;
  }
}

bool SinkCC::SinkCond(int i, int j) {
  const int mask = m_mask_run.as_int(i, j);
  return (mask == 1);
}

void SinkCC::treatInnerMargin(int i, int j,
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

void SinkCC::startNewRun(int i, int j, int &run_number, int &parent, VecList &lists) {
  if (SinkCond(i, j)) {
    parent = 1;
  }
  ConnectedComponents::startNewRun(i, j, run_number, parent, lists);
}

void SinkCC::continueRun(int i, int j, int &run_number, VecList &lists) {
  ConnectedComponents::continueRun(i, j, run_number, lists);
  if (SinkCond(i, j)) {
    setRunSink(run_number, lists["parents"]);
  }
}




FilterExpansionCC::FilterExpansionCC(IceGrid::ConstPtr g, const double fill_value, const IceModelVec2S &bed, const IceModelVec2S &water_level)
  : ValidCC<ConnectedComponents>(g), m_fill_value(fill_value), m_bed(&bed), m_water_level(&water_level) {

  m_min_bed.create(m_grid, "min_bed_mask", WITH_GHOSTS, 1);
  m_min_bed.set(m_fill_value);
  m_max_wl.create(m_grid, "max_water_level", WITH_GHOSTS, 1);
  m_max_wl.set(m_fill_value);

  m_masks.push_back(&m_min_bed);
  m_masks.push_back(&m_max_wl);
  m_fields.push_back(m_bed);
  m_fields.push_back(m_water_level);
}

void FilterExpansionCC::filter_ext(const IceModelVec2S &current_level, const IceModelVec2S &target_level, IceModelVec2Int &mask, IceModelVec2S &min_basin, IceModelVec2S &max_water_level) {
  prepare_mask(current_level, target_level, m_mask_run);
  connected_components::set_validity_mask(4,
                                          m_mask_run,
                                          [](double value) { return value > 1; },
                                          m_mask_validity);

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelMap(run_number, lists, mask, min_basin, max_water_level);
}

void FilterExpansionCC::filter_ext2(const IceModelVec2S &current_level, const IceModelVec2S &target_level, IceModelVec2Int &mask, IceModelVec2S &min_basin, IceModelVec2S &max_water_level) {

  {
    prepare_mask(current_level, target_level, m_mask_run);
    connected_components::set_validity_mask(4,
                                            m_mask_run,
                                            [](double value) { return value > 1; },
                                            m_mask_validity);

    VecList lists;
    unsigned int max_items = 2 * m_grid->ym();
    init_VecList(lists, max_items);

    int run_number = 1;

    compute_runs(run_number, lists, max_items);

    labelMap(run_number, lists, mask, min_basin, max_water_level);
  }

  {
    prepare_mask(target_level, current_level, m_mask_run);
    connected_components::set_validity_mask(4,
                                            m_mask_run,
                                            [](double value) { return value > 1; },
                                            m_mask_validity);

    VecList lists;
    unsigned int max_items = 2 * m_grid->ym();
    init_VecList(lists, max_items);

    int run_number = 1;

    compute_runs(run_number, lists, max_items);

    labelMap2(run_number, lists, mask, min_basin, max_water_level);
  }
}

void FilterExpansionCC::init_VecList(VecList &lists, const unsigned int length) {
  ValidCC<ConnectedComponents>::init_VecList(lists, length);

  std::vector<double> min_bed_list(length);
  std::vector<double> max_wl_list(length);
  lists["min_bed"] = min_bed_list;
  lists["max_wl"]  = max_wl_list;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["min_bed"][k] = m_fill_value;
    lists["max_wl"][k]  = m_fill_value;
  }
}

bool FilterExpansionCC::ForegroundCond(int i, int j) const {
  int mask = m_mask_run.as_int(i, j);
  return mask > 1;
}

void FilterExpansionCC::setRunMinBed(double level, int run, VecList &lists) {
  if (run == 0) {
    return;
  }

  run = connected_components::trackParentRun(run, lists["parents"]);
  if (is_valid(level)) {
    if (is_valid(lists["min_bed"][run])) {
      level = std::min(level, lists["min_bed"][run]);
    }
    lists["min_bed"][run] = level;
  }
}

void FilterExpansionCC::setRunMaxWl(double level, int run, VecList &lists) const {
  if (run == 0) {
    return;
  }

  if (level != m_fill_value) {
    run = connected_components::trackParentRun(run, lists["parents"]);
    if (lists["max_wl"][run] != m_fill_value) {
      level = std::max(level, lists["max_wl"][run]);
    }
    lists["max_wl"][run] = level;
  }
}

void FilterExpansionCC::labelMask(int run_number, const VecList &lists) {
  // set m_mask_run to labels in `lists`:
  connected_components::set_labels(run_number, lists, m_mask_run);

  // use labels in m_mask_run to set other fields:
  IceModelVec::AccessList list{&m_mask_run, &m_mask_validity, &m_min_bed, &m_max_wl};

  const auto
    &validity = lists.find("valid")->second,
    &min_bed   = lists.find("min_bed")->second,
    &max_wl    = lists.find("max_wl")->second;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto label = static_cast<int>(m_mask_run(i, j));

    m_mask_validity(i, j) = validity[label];
    m_min_bed(i, j)       = min_bed[label];
    m_max_wl(i, j)        = max_wl[label];
  }
}

static double min(double a, double b) {
  return (is_valid(a) and ((a < b) or not is_valid(b))) ? a : b;
}

static double max(double a, double b) {
  return (is_valid(a) and ((a > b) or not is_valid(b))) ? a : b;
}

void FilterExpansionCC::treatInnerMargin(int i, int j, bool isNorth, bool isEast, bool isSouth,
                                         bool isWest, VecList &lists, bool &changed) {
  ValidCC<ConnectedComponents>::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  int run = m_mask_run.as_int(i, j);
  if (run <= 0) {
    return;
  }

  auto min_bed = m_min_bed.star(i, j);
  auto max_wl  = m_max_wl.star(i, j);

  double bed_min = min_bed.ij;
  double wl_max  = max_wl.ij;

  StarStencil<bool> flags;
  flags.e = isEast;
  flags.w = isWest;
  flags.n = isNorth;
  flags.s = isSouth;

  for (auto d : { North, East, South, West }) {
    if (flags[d]) {
      bed_min = min(min_bed[d], bed_min);
      wl_max  = max(max_wl[d], wl_max);
    }
  }

  if (bed_min != min_bed.ij) {
    setRunMinBed(bed_min, run, lists);
    changed = true;
  }

  if (wl_max != max_wl.ij) {
    setRunMaxWl(wl_max, run, lists);
    changed = true;
  }
}

void FilterExpansionCC::startNewRun(int i, int j, int &run_number, int &parent, VecList &lists) {
  ValidCC<ConnectedComponents>::startNewRun(i, j, run_number, parent, lists);

  lists["min_bed"][run_number] = (*m_bed)(i, j);
  lists["max_wl"][run_number]  = (*m_water_level)(i, j);
}

void FilterExpansionCC::continueRun(int i, int j, int &run_number, VecList &lists) {
  ValidCC<ConnectedComponents>::continueRun(i, j, run_number, lists);

  setRunMinBed((*m_bed)(i, j), run_number, lists);
  setRunMaxWl((*m_water_level)(i, j), run_number, lists);
}

void FilterExpansionCC::labelMap(int run_number,
                                 const VecList &lists,
                                 IceModelVec2Int &mask,
                                 IceModelVec2S &min_bed,
                                 IceModelVec2S &max_wl) const {

  connected_components::set_labels(run_number, lists, mask);

  const auto
    &valid   = lists.find("valid")->second,
    &bed_min = lists.find("min_bed")->second,
    &wl_max  = lists.find("max_wl")->second;

  IceModelVec::AccessList list{ &mask, &min_bed, &max_wl};

  for (Points p(*mask.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto label = static_cast<int>(mask(i, j));

    if (label > 1) {
      mask(i, j)    = (valid[label] > 0) ? 1 : 2;
      min_bed(i, j) = bed_min[label];
      max_wl(i, j)  = wl_max[label];
    } else {
      mask(i, j)    = 0;
      min_bed(i, j) = m_fill_value;
      max_wl(i, j)  = m_fill_value;
    }
  }
}

void FilterExpansionCC::labelMap2(int run_number, const VecList &lists,
                                  IceModelVec2Int &mask, IceModelVec2S &min_bed,
                                  IceModelVec2S &max_wl) {

  connected_components::set_labels(run_number, lists, mask);

  const auto
    &validity   = lists.find("valid")->second,
    &min_bed_list = lists.find("min_bed")->second,
    &max_wl_list  = lists.find("max_wl")->second;

  IceModelVec::AccessList list{ &mask, &min_bed, &max_wl };

  for (Points p(*mask.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto label = static_cast<int>(mask(i, j));

    if (label > 1) {
      mask(i, j)    = (validity[label] > 0) ? -1 : -2;
      min_bed(i, j) = min_bed_list[label];
      max_wl(i, j)  = max_wl_list[label];
    }
  }
}

/*!
 * Set `result(i, j)` to 2 if `isLake(target_level(i, j))` and not
 * `isLake(current_level(i, j))`, otherwise 0.
 */
void FilterExpansionCC::prepare_mask(const IceModelVec2S &current_level, const IceModelVec2S &target_level,
                                     IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &result, &current_level, &target_level};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (is_valid(target_level(i, j)) and not is_valid(current_level(i, j))) {
      result(i, j) = 2;
    } else {
      result(i, j) = 0;
    }
  }
  result.update_ghosts();
}


namespace connected_components {

/*!
 *
 * Set `result(i, j)` to 1 if `input(i, j) > 1` and at least `threshold` of its neighbors
 * have values greater than 1. Otherwise set `result(i, j)` to 0.
 */
void set_validity_mask(int threshold, const IceModelVec2S &input,
                       const std::function<bool(double)> &condition,
                       IceModelVec2S &result) {

  IceModelVec::AccessList list{ &input, &result };
  for (Points p(*result.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    int n_neighbors = 0;
    if (condition(input(i, j))) {
      auto S = input.star(i, j);
      for (auto direction : { North, East, South, West }) {
        if (condition(S[direction])) {
          ++n_neighbors;
        }
      }
    }

    // Consider a cell to be "valid" if its number of neighbors exceeds the `threshold`.
    if (n_neighbors >= threshold) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  } // end of the loop over grid points

  result.update_ghosts();
}

/*!
 * Replace values in result that satisfy `condition` with `value.`
 */
void replace_values(IceModelVec2S &result,
                    const std::function<bool(double)> &condition, double value) {

  IceModelVec::AccessList list{&result};

  for (Points p(*result.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (condition(result(i, j))) {
      result(i, j) = value;
    }
  }
}

/*!
 * Set the mask in `result` to labels recorded in `lists`, assuming that `lists` has info
 * about `run_number` runs.
 */
void set_labels(int run_number, VecList lists, IceModelVec2S &result) {
  result.set(0.0);

  IceModelVec::AccessList list{&result};

  const auto
    &i_vec   = lists.find("i")->second,
    &j_vec   = lists.find("j")->second,
    &len_vec = lists.find("lengths")->second,
    &parents = lists.find("parents")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    auto j = static_cast<int>(j_vec[k]);
    for (int n = 0; n < len_vec[k]; ++n) {
      auto i = static_cast<int>(i_vec[k]) + n;
      result(i, j) = label;
    }
  }
}

int trackParentRun(int run, const std::vector<double> &parents) {
  while (parents[run] != 0) {
    run = static_cast<int>(parents[run]);
  }
  return run;
}

} // end of namespace connected_components

} //namespace pism
