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

void LakeLevelCC::labelMap(int run_number, const VecList &lists, IceModelVec2S &result) {
  IceModelVec::AccessList list{&result};

  const auto
    &i_vec      = lists.find("i")->second,
    &j_vec      = lists.find("j")->second,
    &len_vec    = lists.find("lengths")->second,
    &parents    = lists.find("parents")->second,
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

bool LakeLevelCC::ForegroundCond(int i, int j) const {
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

  result.set(0);

  labelIsolatedSpots(run_number, lists, result);
}

bool IsolationCC::ForegroundCond(int i, int j) const {
  const double thk = (*m_thk)(i, j);
  const int mask = m_mask_run(i, j);

  return ForegroundCond(thk, mask);
}

void IsolationCC::labelIsolatedSpots(int run_number, const VecList &lists, IceModelVec2Int &result) {
  IceModelVec::AccessList list{&result};
  result.set(0);

  const auto
    &i_vec   = lists.find("i")->second,
    &j_vec   = lists.find("j")->second,
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

void FilterLakesCC::filter_map(int n_filter, IceModelVec2S &lake_level) {
  prepare_mask(lake_level);
  set_mask_validity(n_filter, lake_level);

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelMap(run_number, lists, lake_level);
}

bool FilterLakesCC::ForegroundCond(int i, int j) const {
  const int mask = m_mask_run(i, j);

  return ForegroundCond(mask);
}

void FilterLakesCC::labelMap(int run_number, const VecList &lists, IceModelVec2S &result) {
  IceModelVec::AccessList list{&result};

  const auto
    &i_vec      = lists.find("i")->second,
    &j_vec      = lists.find("j")->second,
    &len_vec    = lists.find("lengths")->second,
    &parents    = lists.find("parents")->second,
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

void FilterLakesCC::set_mask_validity(int n_filter, const IceModelVec2S &lake_level) {
  const Direction dirs[] = { North, East, South, West };

  IceModelVec2S ll_tmp(m_grid, "temp_lake_level", WITH_GHOSTS, 1);
  ll_tmp.copy_from(lake_level);

  IceModelVec::AccessList list{ &m_mask_validity, &ll_tmp };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int n_neighbors = 0;
    if (ll_tmp(i, j) != m_fill_value) {
      StarStencil<double> lake_star = ll_tmp.star(i, j);
      for (int n = 0; n < 4; ++n) {
        const Direction direction = dirs[n];
        if (lake_star[direction] != m_fill_value) {
          ++n_neighbors;
        }
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
                                   const IceModelVec2S &lake_level)
  : ConnectedComponents(g), m_fill_value(fill_value), m_target_level(&target_level),
    m_current_level(&lake_level) {

  m_min_lakelevel.create(m_grid, "min_ll_mask", WITH_GHOSTS, 1);
  m_min_lakelevel.set(m_fill_value);

  m_max_lakelevel.create(m_grid, "max_ll_mask", WITH_GHOSTS, 1);
  m_max_lakelevel.set(m_fill_value);

  m_masks.push_back(&m_min_lakelevel);
  m_masks.push_back(&m_max_lakelevel);

  m_fields.push_back(m_target_level);
  m_fields.push_back(m_current_level);
}

LakePropertiesCC::~LakePropertiesCC() {
  //empty
}

void LakePropertiesCC::getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level) {
  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  min_level.copy_from(m_min_lakelevel);
  max_level.copy_from(m_max_lakelevel);
}

void LakePropertiesCC::init_VecList(VecList &lists, const unsigned int length) {
  ConnectedComponents::init_VecList(lists, length);

  std::vector<double> min_ll_list(length), max_ll_list(length);
  lists["min_ll"]  = min_ll_list;
  lists["max_ll"]  = max_ll_list;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["min_ll"][k]  = m_fill_value;
    lists["max_ll"][k]  = m_fill_value;
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


bool LakePropertiesCC::ForegroundCond(int i, int j) const {
  const double target  = (*m_target_level)(i, j);

  return isLake(target);
}

void LakePropertiesCC::labelMask(int run_number, const VecList &lists) {
  IceModelVec::AccessList list;
  addFieldVecAccessList(m_masks, list);

  const auto
    &i_vec   = lists.find("i")->second,
    &j_vec   = lists.find("j")->second,
    &len_vec = lists.find("lengths")->second,
    &parents = lists.find("parents")->second,
    &min_ll  = lists.find("min_ll")->second,
    &max_ll  = lists.find("max_ll")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    const double min_ll_label = min_ll[label],
                 max_ll_label = max_ll[label];
    unsigned int j = j_vec[k];
    for (unsigned int n = 0; n < len_vec[k]; ++n) {
      const int i = i_vec[k] + n;

      m_mask_run(i, j) = label;
      m_min_lakelevel(i, j) = min_ll_label;
      m_max_lakelevel(i, j) = max_ll_label;
    }
  }
}

void LakePropertiesCC::treatInnerMargin(int i, int j,
                                        const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                        VecList &lists, bool &changed) {
  ConnectedComponents::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  int run = m_mask_run.as_int(i, j);
  if (run > 0) {
    StarStencil<double> min_ll_star  = m_min_lakelevel.star(i, j),
                        max_ll_star  = m_max_lakelevel.star(i, j);

    double min_level = min_ll_star.ij,
           max_level = max_ll_star.ij;

    if (isWest) {
      if (isLake(min_ll_star.w) and ((min_ll_star.w < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.w;
      }
      if (isLake(max_ll_star.w) and ((max_ll_star.w > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.w;
      }
    }
    if (isNorth) {
      if (isLake(min_ll_star.n) and ((min_ll_star.n < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.n;
      }
      if (isLake(max_ll_star.n) and ((max_ll_star.n > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.n;
      }
    }
    if (isEast) {
      if (isLake(min_ll_star.e) and ((min_ll_star.e < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.e;
      }
      if (isLake(max_ll_star.e) and ((max_ll_star.e > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.e;
      }
    }
    if (isSouth) {
      if (isLake(min_ll_star.s) and ((min_ll_star.s < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.s;
      }
      if (isLake(max_ll_star.s) and ((max_ll_star.s > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.s;
      }
    }
    if (min_level != min_ll_star.ij) {
      setRunMinLevel(min_level, run, lists);
      changed = true;
    }
    if (max_level != max_ll_star.ij) {
      setRunMaxLevel(max_level, run, lists);
      changed = true;
    }
  }
}

void LakePropertiesCC::startNewRun(int i, int j, int &run_number, int &parent, VecList &lists) {
  ConnectedComponents::startNewRun(i, j, run_number, parent, lists);

  lists["min_ll"][run_number]  = (*m_current_level)(i, j);
  lists["max_ll"][run_number]  = (*m_current_level)(i, j);
}

void LakePropertiesCC::continueRun(int i, int j, int &run_number, VecList &lists) {
  ConnectedComponents::continueRun(i, j, run_number, lists);

  setRunMinLevel((*m_current_level)(i, j), run_number, lists);
  setRunMaxLevel((*m_current_level)(i, j), run_number, lists);
}


LakeAccumulatorCCSerial::LakeAccumulatorCCSerial(IceGrid::ConstPtr g, const double fill_value)
  : ConnectedComponentsSerial(g),
    m_fill_value(fill_value),
    m_initialized(false) {
  //empty
}

LakeAccumulatorCCSerial::~LakeAccumulatorCCSerial() {
  //empty
}

void LakeAccumulatorCCSerial::init(const IceModelVec2S &lake_level) {
  prepare_mask(lake_level);

  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(m_lists, max_items);

  m_run_number = 0;

  compute_runs(m_run_number, m_lists, max_items);

  m_initialized = true;
}

void LakeAccumulatorCCSerial::accumulate(const IceModelVec2S &in, IceModelVec2S &result) {

  if (not m_initialized) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "LakeAccumulatorCCSerial is not initialized.");
  }

  petsc::Vec::Ptr in_vec_p0 = in.allocate_proc0_copy(),
                  result_vec_p0 = result.allocate_proc0_copy();
  in.put_on_proc0(*in_vec_p0);

  //Init result and put it on proc0
  result.set(m_fill_value);
  result.put_on_proc0(*result_vec_p0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      petsc::VecArray2D in_p0(*in_vec_p0, m_grid->Mx(), m_grid->My()),
                        result_p0(*result_vec_p0, m_grid->Mx(), m_grid->My());
      //Init allocator
      std::vector<double> accumulator(m_run_number + 1, 0.0);

      const auto
        &i_vec   = m_lists.find("i")->second,
        &j_vec   = m_lists.find("j")->second,
        &len_vec = m_lists.find("lengths")->second,
        &parents = m_lists.find("parents")->second;

      //accumulate values
      for (int k = 0; k <= m_run_number; ++k) {
        const int j = j_vec[k];
        const int label = trackParentRun(k, parents);
        for (unsigned int n = 0; n < len_vec[k]; ++n) {
          const int i = i_vec[k] + n;
          accumulator[label] += in_p0(i, j);
        }
      }

      //label result
      for (int k = 0; k <= m_run_number; ++k) {
        const int j = j_vec[k];
        const int label = trackParentRun(k, parents);
        for (unsigned int n = 0; n < len_vec[k]; ++n) {
          const int i = i_vec[k] + n;
          result_p0(i, j) = accumulator[label];
        }
      }
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  //Get result from Processor 0
  result.get_from_proc0(*result_vec_p0);
}

bool LakeAccumulatorCCSerial::ForegroundCond(int i, int j) const {
  const int mask = (*m_mask_run_p0_ptr)(i, j);
  return ForegroundCond(mask);
}

void LakeAccumulatorCCSerial::prepare_mask(const IceModelVec2S &lake_level) {
  IceModelVec::AccessList list{ &m_mask_run, &lake_level};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (isLake(lake_level(i, j))) {
      m_mask_run(i, j) = 1;
    } else {
      m_mask_run(i, j) = 0;
    }
  }
}

} // namespace pism
