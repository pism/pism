#include "SeaLevel_ConnectedComponents.hh"

namespace pism{

SeaLevelCC::SeaLevelCC(IceGrid::ConstPtr g,
                       const double drho,
                       const IceModelVec2S &bed,
                       const IceModelVec2S &thk,
                       const double fill_value)
  :FillingAlgCC<SinkCC>(g, drho, bed, thk, fill_value) {

  // prepare the mask
  {
    IceModelVec::AccessList list{ &m_mask_run };
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Set "sink" at a margin of the computational domain
      m_mask_run(i, j) = grid_edge(*m_grid, i, j) ? 1 : 0;
    }
    m_mask_run.update_ghosts();
  }
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

  // Initialize the mask:
  result.set(0.0);
  connected_components::set_labels(run_number, lists, result);
  connected_components::replace_values(result,
                                       [](double label) { return label > 1; },
                                       1);
}

bool SeaLevelCC::ForegroundCond(int i, int j) const {
  double bed = (*m_bed)(i, j),
         thk = (*m_thk)(i, j),
         sea_level = (*m_sea_level)(i, j);
  int mask = m_mask_run.as_int(i, j);

  return is_foreground(bed, thk, mask, sea_level, m_offset);
}

}
