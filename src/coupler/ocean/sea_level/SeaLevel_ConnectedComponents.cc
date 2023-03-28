#include "SeaLevel_ConnectedComponents.hh"

namespace pism{

SeaLevelCC::SeaLevelCC(IceGrid::ConstPtr g,
                       const double drho,
                       const IceModelVec2S &bed,
                       const IceModelVec2S &thk,
                       const double fill_value)
  :FillingAlgCC<SinkCC>(g, drho, bed, thk, fill_value) {
  prepare_mask();
}

SeaLevelCC::SeaLevelCC(IceGrid::ConstPtr g,
                       const double drho,
                       const IceModelVec2S &bed,
                       const IceModelVec2S &thk,
                       const IceModelVec2Int &run_mask,
                       const double fill_value)
  :FillingAlgCC<SinkCC>(g, drho, bed, thk, fill_value) {
  m_mask_run.copy_from(run_mask);
}

SeaLevelCC::~SeaLevelCC() {
  //empty
}

void SeaLevelCC::computeSeaLevel(IceModelVec2S &SeaLevel, const double Offset) {
  computeSeaLevel(SeaLevel, Offset, SeaLevel);
}

void SeaLevelCC::computeSeaLevel(const IceModelVec2S &SeaLevel, const double Offset, IceModelVec2S &result) {
  m_sea_level = &SeaLevel;
  m_offset = Offset;

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  IceModelVec::AccessList list{ m_sea_level };
  compute_runs(run_number, lists, max_items);

  labelSLMap(run_number, lists, result);
}

void SeaLevelCC::computeMask(const IceModelVec2S &SeaLevel, const double Offset, IceModelVec2Int &result) {
  m_sea_level = &SeaLevel;
  m_offset = Offset;

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  IceModelVec::AccessList list{ m_sea_level };
  compute_runs(run_number, lists, max_items);

  labelSLMask(run_number, lists, result);
}

void SeaLevelCC::labelSLMap(int run_number, const VecList lists, IceModelVec2S &result) {
  IceModelVec::AccessList list{&result};

  const auto
    &i_vec   = lists.find("i")->second,
    &j_vec   = lists.find("j")->second,
    &len_vec = lists.find("lengths")->second,
    &parents = lists.find("parents")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    if(label > 1) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = m_fill_value;
      }
    }
  }
}

void SeaLevelCC::labelSLMask(int run_number, const VecList lists, IceModelVec2Int &result) {
  //Init mask to 0
  result.set(0);

  IceModelVec::AccessList list{&result};

  const auto
    &i_vec   = lists.find("i")->second,
    &j_vec   = lists.find("j")->second,
    &len_vec = lists.find("lengths")->second,
    &parents = lists.find("parents")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    if(label > 1) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = 1;
      }
    }
  }
}

void SeaLevelCC::prepare_mask() {
  IceModelVec::AccessList list{&m_mask_run};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Set "sink" at a margin of the computational domain
    m_mask_run(i, j) = grid_edge(*m_grid, i, j) ? 1 : 0;
  }
  m_mask_run.update_ghosts();
}

bool SeaLevelCC::ForegroundCond(int i, int j) const {
  double bed = (*m_bed)(i, j),
         thk = (*m_thk)(i, j),
         sea_level = (*m_sea_level)(i, j);
  int mask = m_mask_run(i, j);

  return is_foreground(bed, thk, mask, sea_level, m_offset);
}

}
