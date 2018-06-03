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
               &len_vec    = lists.find("lengths")->second,
               &parents    = lists.find("parents")->second;

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

} // namespace pism
