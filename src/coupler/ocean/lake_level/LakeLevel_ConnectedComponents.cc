#include "LakeLevel_ConnectedComponents.hh"

namespace pism {

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, const double fill_value)
  : FillingAlgCC<ValidCC<SinkCC> >(g, drho, bed, thk, fill_value) {
  IceModelVec2CellType pism_mask_type;
  pism_mask_type.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  pism_mask_type.copy_from(pism_mask);
  prepare_mask(pism_mask_type);
  m_mask_validity.set(1);
}

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, const double fill_value,
                         const IceModelVec2Int &valid_mask)
  : FillingAlgCC<ValidCC<SinkCC> >(g, drho, bed, thk, fill_value) {
  IceModelVec2CellType pism_mask_type;
  pism_mask_type.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  pism_mask_type.copy_from(pism_mask);
  this->prepare_mask(pism_mask_type);
  m_mask_validity.copy_from(valid_mask);
}

LakeLevelCC::~LakeLevelCC() {
  //empty
}

void LakeLevelCC::computeLakeLevel(const double zMin, const double zMax, const double dz, const double offset, IceModelVec2S &result) {
  m_offset = offset;
  result.set(m_fill_value);

  double lakeLevel = zMin;
  while (lakeLevel <= zMax) {
    fill2Level(lakeLevel, result);
    lakeLevel += dz;
  }
}

void LakeLevelCC::fill2Level(const double level, IceModelVec2S &result) {
  m_level = level;

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelMap(run_number, lists, result);
}

void LakeLevelCC::labelMap(const int run_number, const VecList &lists, IceModelVec2S &result) {
  IceModelVec::AccessList list{&result};

  const RunVec &i_vec = lists.find("i")->second,
               &j_vec = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second,
               &valid_list = lists.find("valid")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    const bool validLake = ((label > 1) and (valid_list[label] > 0));
    if (validLake) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = m_level;
      }
    }
  }
}

void LakeLevelCC::prepare_mask(const IceModelVec2CellType &pism_mask) {
  IceModelVec::AccessList list{ &m_mask_run, &pism_mask };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest = (i == m_i_global_first), isEast = (i == m_i_global_last), isSouth = (j == m_j_global_first),
         isNorth = (j == m_j_global_last), isMargin = (isWest or isEast or isSouth or isNorth);

    //Set sink, where pism_mask is ocean or at margin of computational domain
    if (isMargin or pism_mask.ocean(i, j)) {
      m_mask_run(i, j) = 1;
    } else {
      m_mask_run(i, j) = 0;
    }
  }
  m_mask_run.update_ghosts();
}

bool LakeLevelCC::ForegroundCond(const int i, const int j) const {
  double bed = (*m_bed)(i, j),
         thk = (*m_thk)(i, j);
  int mask = m_mask_run(i, j);

  return FillingAlgCC::ForegroundCond(bed, thk, mask, m_level, m_offset);
}


IsolationCC::IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk,
                         const double thk_theshold)
  : SinkCC(g), m_thk_threshold(thk_theshold), m_thk(&thk) {
  prepare_mask();
  m_fields.push_back(m_thk);
}

IsolationCC::~IsolationCC() {
  //empty
}

void IsolationCC::find_isolated_spots(IceModelVec2Int &result) {
  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelIsolatedSpots(run_number, lists, result);
}

bool IsolationCC::ForegroundCond(const int i, const int j) const {
  const double thk = (*m_thk)(i, j);
  const int mask = m_mask_run(i, j);

  return ForegroundCond(thk, mask);
}

void IsolationCC::labelIsolatedSpots(const int run_number, const VecList &lists, IceModelVec2Int &result) {
  IceModelVec::AccessList list{&result};
  result.set(0);

  const RunVec &i_vec = lists.find("i")->second,
               &j_vec = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    if (label == 1) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = 1;
      }
    }
  }
}

void IsolationCC::prepare_mask() {
  IceModelVec::AccessList list{ &m_mask_run };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest = (i == m_i_global_first), isEast = (i == m_i_global_last), isSouth = (j == m_j_global_first),
         isNorth = (j == m_j_global_last), isMargin = (isWest or isEast or isSouth or isNorth);

    //Set not isolated at margin
    m_mask_run(i, j) = isMargin ? 1 : 0;
  }
  m_mask_run.update_ghosts();
}


FilterLakesCC::FilterLakesCC(IceGrid::ConstPtr g, const double fill_value)
  : ValidCC<ConnectedComponents>(g), m_fill_value(fill_value) {
  //empty
}

FilterLakesCC::~FilterLakesCC() {

}

void FilterLakesCC::filter_map(const int n_filter, IceModelVec2S &lake_level) {
  prepare_mask(lake_level);
  set_mask_validity(n_filter, lake_level);

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelMap(run_number, lists, lake_level);
}

bool FilterLakesCC::ForegroundCond(const int i, const int j) const {
  const int mask = m_mask_run(i, j);

  return ForegroundCond(mask);
}

void FilterLakesCC::labelMap(const int run_number, const VecList &lists, IceModelVec2S &result) {
  IceModelVec::AccessList list{&result};

  const RunVec &i_vec = lists.find("i")->second,
               &j_vec = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second,
               &valid_list = lists.find("valid")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    const bool valid = (valid_list[label] > 0);
    if (not valid) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = m_fill_value;
      }
    }
  }
}

void FilterLakesCC::prepare_mask(const IceModelVec2S &lake_level) {
  IceModelVec::AccessList list{ &m_mask_run, &lake_level};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    //Set sink, where pism_mask is ocean or at margin of computational domain
    if (isLake(lake_level(i, j))) {
      m_mask_run(i, j) = 2;
    } else {
      m_mask_run(i, j) = 0;
    }
  }
  m_mask_run.update_ghosts();
}

void FilterLakesCC::set_mask_validity(const int n_filter, const IceModelVec2S &lake_level) {
  const Direction dirs[] = { North, East, South, West };

  IceModelVec2S ll_tmp(m_grid, "temp_lake_level", WITH_GHOSTS, 1);
  ll_tmp.copy_from(lake_level);

  IceModelVec::AccessList list{ &m_mask_validity, &ll_tmp };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    StarStencil<double> lake_star = ll_tmp.star(i, j);

    int n_neighbors = 0;
    for (int n = 0; n < 4; ++n) {
      const Direction direction = dirs[n];
      if (lake_star[direction] != m_fill_value) {
        ++n_neighbors;
      }
    }
    //Set cell valid if number of neighbors exceeds threshold
    if (n_neighbors >= n_filter) {
      m_mask_validity(i, j) = 1;
    } else {
      m_mask_validity(i, j) = 0;
    }
  }
  m_mask_validity.update_ghosts();
}


LakePropertiesCC::LakePropertiesCC(IceGrid::ConstPtr g, const double fill_value, const IceModelVec2S &target_level,
                                   const IceModelVec2S &lake_level, const IceModelVec2S &floating_thresh)
  : ConnectedComponents(g), m_fill_value(fill_value), m_target_level(&target_level),
    m_current_level(&lake_level), m_floating_threshold_level(&floating_thresh) {

  m_min_ll_tmp.create(m_grid, "min_ll_mask", WITH_GHOSTS, 1);
  m_min_ll_tmp.set(m_fill_value);

  m_max_ll_tmp.create(m_grid, "max_ll_mask", WITH_GHOSTS, 1);
  m_max_ll_tmp.set(m_fill_value);

  m_min_float_tmp.create(m_grid, "min_float_mask", WITH_GHOSTS, 1);
  m_min_float_tmp.set(m_fill_value);

  m_masks.push_back(&m_min_ll_tmp);
  m_masks.push_back(&m_max_ll_tmp);
  m_masks.push_back(&m_min_float_tmp);

  m_fields.push_back(m_target_level);
  m_fields.push_back(m_current_level);
  m_fields.push_back(m_floating_threshold_level);
}

LakePropertiesCC::~LakePropertiesCC() {
  //empty
}

void LakePropertiesCC::getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level,
                                         IceModelVec2S &min_float_level) {
  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  min_level.copy_from(m_min_ll_tmp);
  max_level.copy_from(m_max_ll_tmp);
  min_float_level.copy_from(m_min_float_tmp);
}

void LakePropertiesCC::init_VecList(VecList &lists, const unsigned int length) {
  ConnectedComponents::init_VecList(lists, length);

  RunVec min_ll_list(length), max_ll_list(length), min_float_list(length);
  lists["min_ll"] = min_ll_list;
  lists["max_ll"] = max_ll_list;
  lists["min_float"] = min_float_list;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["min_ll"][k] = m_fill_value;
    lists["max_ll"][k] = m_fill_value;
    lists["min_float"][k] = m_fill_value;
  }
}

void LakePropertiesCC::setRunMinLevel(double level, int run, VecList &lists) {
  if (run == 0) {
    return;
  }

  run = trackParentRun(run, lists["parents"]);
  if (isLake(level)) {
    if (isLake(lists["min_ll"][run])) {
      level = std::min(level, lists["min_ll"][run]);
    }
    lists["min_ll"][run] = level;
  }
}

void LakePropertiesCC::setRunMaxLevel(double level, int run, VecList &lists) {
  if (run == 0) {
    return;
  }

  run = trackParentRun(run, lists["parents"]);
  if (isLake(level)) {
    if (isLake(lists["max_ll"][run])) {
      level = std::max(level, lists["max_ll"][run]);
    }
    lists["max_ll"][run] = level;
  }
}

void LakePropertiesCC::setRunMinFloatLevel(double level, int run, VecList &lists) {
  if (run == 0) {
    return;
  }

  run = trackParentRun(run, lists["parents"]);
  if (isLake(level)) {
    if (isLake(lists["min_float"][run])) {
      level = std::min(level, lists["min_float"][run]);
    }
    lists["min_float"][run] = level;
  }
}

bool LakePropertiesCC::ForegroundCond(const int i, const int j) const {
  const double target  = (*m_target_level)(i, j),
               current = (*m_current_level)(i, j);

  return (ForegroundCond(target, current));
}

void LakePropertiesCC::labelMask(int run_number, const VecList &lists) {
  IceModelVec::AccessList list;
  addFieldVecAccessList(m_masks, list);

  const RunVec &i_vec = lists.find("i")->second,
               &j_vec = lists.find("j")->second,
               &len_vec = lists.find("lengths")->second,
               &parents = lists.find("parents")->second,
               &min_vec = lists.find("min_ll")->second,
               &max_vec = lists.find("max_ll")->second,
               &min_float_vec = lists.find("min_float")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    const double min_ll_label = min_vec[label],
                 max_ll_label = max_vec[label],
                 min_float_label = min_float_vec[label];
    unsigned int j = j_vec[k];
    for (unsigned int n = 0; n < len_vec[k]; ++n) {
      const int i = i_vec[k] + n;

      m_mask_run(i, j) = label;
      m_min_ll_tmp(i, j) = min_ll_label;
      m_max_ll_tmp(i, j) = max_ll_label;
      m_min_float_tmp(i, j) = min_float_label;
    }
  }
}

void LakePropertiesCC::treatInnerMargin(const int i, const int j,
                                        const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                        VecList &lists, bool &changed) {
  ConnectedComponents::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  int run = m_mask_run.as_int(i, j);
  if (run > 0) {
    StarStencil<double> min_star = m_min_ll_tmp.star(i, j),
                        max_star = m_max_ll_tmp.star(i, j),
                        minfloat_star = m_min_float_tmp.star(i, j);

    double min_level = min_star.ij,
           max_level = max_star.ij,
           min_float = minfloat_star.ij;

    if (isWest) {
      if (isLake(min_star.w) and ((min_star.w < min_level) or not isLake(min_level))) {
        min_level = min_star.w;
      }
      if (isLake(max_star.w) and ((max_star.w > max_level) or not isLake(max_level))) {
        max_level = max_star.w;
      }
      if (isLake(minfloat_star.w) and ((minfloat_star.w < min_float) or not isLake(min_float))) {
        min_float = minfloat_star.w;
      }
    }
    if (isNorth) {
      if (isLake(min_star.n) and ((min_star.n < min_level) or not isLake(min_level))) {
        min_level = min_star.n;
      }
      if (isLake(max_star.n) and ((max_star.n > max_level) or not isLake(max_level))) {
        max_level = max_star.n;
      }
      if (isLake(minfloat_star.n) and ((minfloat_star.n < min_float) or not isLake(min_float))) {
        min_float = minfloat_star.n;
      }
    }
    if (isEast) {
      if (isLake(min_star.e) and ((min_star.e < min_level) or not isLake(min_level))) {
        min_level = min_star.e;
      }
      if (isLake(max_star.e) and ((max_star.e > max_level) or not isLake(max_level))) {
        max_level = max_star.e;
      }
      if (isLake(minfloat_star.e) and ((minfloat_star.e < min_float) or not isLake(min_float))) {
        min_float = minfloat_star.e;
      }
    }
    if (isSouth) {
      if (isLake(min_star.s) and ((min_star.s < min_level) or not isLake(min_level))) {
        min_level = min_star.s;
      }
      if (isLake(max_star.s) and ((max_star.s > max_level) or not isLake(max_level))) {
        max_level = max_star.s;
      }
      if (isLake(minfloat_star.s) and ((minfloat_star.s < min_float) or not isLake(min_float))) {
        min_float = minfloat_star.s;
      }
    }
    if (min_level != min_star.ij) {
      setRunMinLevel(min_level, run, lists);
      changed = true;
    }
    if (max_level != max_star.ij) {
      setRunMaxLevel(max_level, run, lists);
      changed = true;
    }
    if (min_float != minfloat_star.ij) {
      setRunMinFloatLevel(min_float, run, lists);
      changed = true;
    }
  }
}

void LakePropertiesCC::startNewRun(const int i, const int j, int &run_number, int &parent, VecList &lists) {
  ConnectedComponents::startNewRun(i, j, run_number, parent, lists);

  lists["min_ll"][run_number] = (*m_current_level)(i, j);
  lists["max_ll"][run_number] = (*m_current_level)(i, j);
  lists["min_float"][run_number] = (*m_floating_threshold_level)(i, j);
}

void LakePropertiesCC::continueRun(const int i, const int j, int &run_number, VecList &lists) {
  ConnectedComponents::continueRun(i, j, run_number, lists);

  setRunMinLevel((*m_current_level)(i, j), run_number, lists);
  setRunMaxLevel((*m_current_level)(i, j), run_number, lists);
  setRunMinFloatLevel((*m_floating_threshold_level)(i, j), run_number, lists);
}

} // namespace pism