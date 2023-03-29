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
    // Set "sink" if pism_mask is ocean or at a margin of the computational domain
    if (grid_edge(*m_grid, i, j) or pism_mask.ocean(i, j)) {
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
  int mask = m_mask_run.as_int(i, j);

  return is_foreground(bed, thk, mask, m_level, m_offset);
}



LakeAccumulatorCCSerial::LakeAccumulatorCCSerial(IceGrid::ConstPtr g, const double fill_value)
  : ConnectedComponentsSerial(g),
    m_initialized(false),
    m_fill_value(fill_value)
{
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
  return (mask > 0);
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
