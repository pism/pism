#include "LakeLevel_ConnectedComponents.hh"

namespace pism {

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, const double drho,
                         const double thk_threshold, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
                         const double fill_value)
    : FillingAlgCC(g, drho, bed, thk, true, fill_value),
      m_thk_threshold(thk_threshold) {

  this->init(pism_mask);
  this->prepare_mask();
  m_mask_validity.set(1);
}

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, const double drho,
                         const double thk_threshold, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
                         IceModelVec2Int &run_mask, const double fill_value)
    : FillingAlgCC(g, drho, bed, thk, true, fill_value),
      m_thk_threshold(thk_threshold) {

  this->init(pism_mask);
  m_mask_run.copy_from(run_mask);
  m_mask_validity.set(1);
}

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, const double drho,
                         const double thk_threshold, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
                         const double fill_value, IceModelVec2Int &valid_mask)
    : FillingAlgCC(g, drho, bed, thk, true, fill_value),
      m_thk_threshold(thk_threshold) {

  this->init(pism_mask);
  this->prepare_mask();
  m_mask_validity.copy_from(valid_mask);
}

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, const double drho,
                         const double thk_threshold, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, IceModelVec2Int &pism_mask,
                         IceModelVec2Int &run_mask, const double fill_value,
                         IceModelVec2Int &valid_mask)
    : FillingAlgCC(g, drho, bed, thk, true, fill_value),
      m_thk_threshold(thk_threshold) {

  this->init(pism_mask);
  m_mask_run.copy_from(run_mask);
  m_mask_validity.copy_from(valid_mask);
}

LakeLevelCC::~LakeLevelCC() {
  //empty
}

void LakeLevelCC::init(IceModelVec2Int &pism_mask) {
  m_pism_mask.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  m_pism_mask.copy_from(pism_mask);
  m_floatation_level.set(m_fill_value);
}

void LakeLevelCC::floodMap(double zMin, double zMax, double dz) {
  double lakeLevel = zMin;
  while (lakeLevel <= zMax) {
    this->fill2Level(lakeLevel);
    lakeLevel += dz;
  }
}

void LakeLevelCC::labelMap_impl(unsigned int run_number, std::vector<unsigned int> &i_vec,
                                std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents,
                                std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) {
  IceModelVec::AccessList list{ &m_floatation_level, &m_level };
  for (unsigned int k = 0; k <= run_number; ++k) {
    unsigned int label = trackParentRun(k, parents);
    const bool isLake_label = ((label > 1) and isValidList[label]);
    for (unsigned int n = 0; n < lengths[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      if (isLake_label) {
        m_floatation_level(i, j) = m_level(i, j);
      }
    }
  }
}

void LakeLevelCC::lake_levels(IceModelVec2S &result) const {
  this->get_floatation_level(result);
}

void LakeLevelCC::prepare_mask_impl() {

  IceModelVec::AccessList list{ &m_mask_run, &m_pism_mask };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest = (i == m_i_global_first), isEast = (i == m_i_global_last), isSouth = (j == m_j_global_first),
         isNorth = (j == m_j_global_last), isMargin = (isWest or isEast or isSouth or isNorth);

    //Set sink, where pism_mask is ocean or at margin of computational domain
    if (isMargin or m_pism_mask.ocean(i, j)) {
      m_mask_run(i, j) = 1;
    } else {
      m_mask_run(i, j) = 0;
    }
  }
  m_mask_run.update_ghosts();
}



IsolationCC::IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk, const double thk_threshold)
  :FillingAlgCC(g, 0.0, thk, thk, false, 0.0), m_thk_threshold(thk_threshold){
  this->prepare_mask();
}

IsolationCC::~IsolationCC() {
  //empty
}

void IsolationCC::find_isolated_spots() {
  this->fill2Level(0.0);
}

void IsolationCC::isolation_mask(IceModelVec2Int &result) const {
  this->get_floatation_level(result);
}

bool IsolationCC::ForegroundCond_impl(double bed, double thk, int mask, double Level, double Offset) const {
  (void) bed;
  (void) Level;
  (void) Offset;
  return ((thk < m_thk_threshold) or (mask > 0));
}

void IsolationCC::labelMap_impl(unsigned int run_number, std::vector<unsigned int> &i_vec, std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) {
  (void) isValidList;
  //Set mask to 1, where cell is not isolated by ice
  IceModelVec::AccessList list{&m_floatation_level};
  m_floatation_level.set(0.0);

  for(unsigned int k = 0; k <= run_number; ++k) {
    unsigned int label = trackParentRun(k, parents);
    if(label == 1) {
      for(unsigned int n = 0; n < lengths[k]; ++n) {
        const int i = i_vec[k] + n,
                  j = j_vec[k];
        m_floatation_level(i, j) = 1;
      }
    }
  }
}

void IsolationCC::prepare_mask_impl() {

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



FilterLakesCC::FilterLakesCC(IceGrid::ConstPtr g, const IceModelVec2S &lake_levels, const double fill_value)
  :FillingAlgCC(g, 0.0, lake_levels, lake_levels, true, fill_value) {
  this->prepare_mask();
  //The FillingAlgCC class is modified to be used in the filter algorithm
  //In FillingAlgCC the variable m_bed holds the  lake level variable
  //To stop confusion m_lake_level is used here and points to m_bed
  m_lake_level = m_bed;
  prepare_mask();
}

FilterLakesCC::~FilterLakesCC() {
  //empty
}

void FilterLakesCC::set_mask_validity(int n_filter) {
  const Direction dirs[] = { North, East, South, West };

  IceModelVec2S sl_tmp(m_grid, "temp_lake_level", WITH_GHOSTS, 1);
  sl_tmp.copy_from(*m_lake_level);

  IceModelVec::AccessList list{ &m_mask_validity, &sl_tmp };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    StarStencil<double> lake_star = sl_tmp.star(i, j);

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

void FilterLakesCC::filter_map(int n_filter) {
  set_mask_validity(n_filter);
  this->fill2Level();
}

void FilterLakesCC::filtered_levels(IceModelVec2S &result) const {
  this->get_floatation_level(result);
}

bool FilterLakesCC::ForegroundCond_impl(double lake_level, double thk, int mask, double Level, double Offset) const {
  (void) thk;
  (void) mask;
  (void) Level;
  (void) Offset;
  return (lake_level != m_fill_value);
}

void FilterLakesCC::labelMap_impl(unsigned int run_number, std::vector<unsigned int> &i_vec, std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) {
  //Set mask to 1, where cell is not isolated by ice
  IceModelVec::AccessList list{&m_floatation_level, m_lake_level};
  m_floatation_level.set(m_fill_value);

  for(unsigned int k = 0; k <= run_number; ++k) {
    unsigned int label = trackParentRun(k, parents);
    bool isValid = isValidList[label];
    if(isValid) {
      for(unsigned int n = 0; n < lengths[k]; ++n) {
        const int i = i_vec[k] + n,
                  j = j_vec[k];
        m_floatation_level(i, j) = (*m_lake_level)(i, j);
      }
    }
  }
}

void FilterLakesCC::prepare_mask_impl() {
  m_mask_run.set(0);
}



} // namespace pism
