#include "LakeLevel_ConnectedComponents.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, double drho, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, const IceModelVec2Int &pism_mask, double fill_value)
  : FillingAlgCC<ValidCC<SinkCC> >(g, drho, bed, thk, fill_value) {
  IceModelVec2CellType pism_mask_type;
  pism_mask_type.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  pism_mask_type.copy_from(pism_mask);
  prepare_mask(pism_mask_type, m_mask_run);
  m_mask_validity.set(1);
}

LakeLevelCC::LakeLevelCC(IceGrid::ConstPtr g, double drho, const IceModelVec2S &bed,
                         const IceModelVec2S &thk, const IceModelVec2Int &pism_mask,
                         double fill_value, const IceModelVec2Int &valid_mask)
  : FillingAlgCC<ValidCC<SinkCC> >(g, drho, bed, thk, fill_value) {
  IceModelVec2CellType pism_mask_type;
  pism_mask_type.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  pism_mask_type.copy_from(pism_mask);
  this->prepare_mask(pism_mask_type, m_mask_run);
  m_mask_validity.copy_from(valid_mask);
}

void LakeLevelCC::computeLakeLevel(double zMin, double zMax, double dz, double offset, IceModelVec2S &result) {
  m_offset = offset;
  result.set(m_fill_value);

  double z = zMin;
  while (z <= zMax) {
    fill2Level(z, result);
    z += dz;
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

void LakeLevelCC::labelMap(int run_number, const VecList &lists, IceModelVec2S &result) const {
  IceModelVec::AccessList list{&result};

  const auto
    &i_vec      = lists.find("i")->second,
    &j_vec      = lists.find("j")->second,
    &len_vec    = lists.find("lengths")->second,
    &parents    = lists.find("parents")->second,
    &valid_list = lists.find("valid")->second;

  for (int k = 0; k <= run_number; ++k) {
    int label      = connected_components::trackParentRun(k, parents);
    bool validLake = ((label > 1) and (valid_list[label] > 0));
    if (validLake) {
      const int j = j_vec[k];
      for (int n = 0; n < len_vec[k]; ++n) {
        const int i  = i_vec[k] + n;
        result(i, j) = m_level;
      }
    }
  }
}

void LakeLevelCC::prepare_mask(const IceModelVec2CellType &pism_mask, IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &result, &pism_mask };
  for (Points p(*m_grid); p; p.next()) {

    const int i = p.i(), j = p.j();
    // Set "sink" if pism_mask is ocean or at a margin of the computational domain
    if (grid_edge(*m_grid, i, j) or pism_mask.ocean(i, j)) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  }
  result.update_ghosts();
}

bool LakeLevelCC::ForegroundCond(int i, int j) const {
  return is_foreground((*m_bed)(i, j), (*m_thk)(i, j), m_mask_run.as_int(i, j), m_level, m_offset);
}



} // namespace pism
