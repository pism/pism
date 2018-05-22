#include "SeaLevel_ConnectedComponents.hh"

namespace pism{

SeaLevelCC::SeaLevelCC(IceGrid::ConstPtr g,
                       const double drho,
                       const IceModelVec2S &bed,
                       const IceModelVec2S &thk,
                       const double fill_value)
  :FillingAlgCC(g, drho, bed, thk, false, fill_value) {
  this->prepare_mask();
  m_crop_mask.create(m_grid, "sl_crop_mask", WITHOUT_GHOSTS);
  m_crop_mask.set(0);
}

SeaLevelCC::SeaLevelCC(IceGrid::ConstPtr g,
                       const double drho,
                       const IceModelVec2S &bed,
                       const IceModelVec2S &thk,
                       const IceModelVec2Int &run_mask,
                       const double fill_value)
  :FillingAlgCC(g, drho, bed, thk, false, fill_value) {
  m_mask_run.copy_from(run_mask);
  m_crop_mask.create(m_grid, "sl_crop_mask", WITHOUT_GHOSTS);
  m_crop_mask.set(0);
}

SeaLevelCC::~SeaLevelCC() {
  //empty
}

void SeaLevelCC::fill2SeaLevel(double SeaLevel) {
  m_floatation_level.set(SeaLevel);
  this->fill2Level(SeaLevel);
}

void SeaLevelCC::fill2SeaLevel(double SeaLevel, double Offset) {
  m_floatation_level.set(SeaLevel);
  this->fill2Level(SeaLevel, Offset);
}

void SeaLevelCC::fill2SeaLevel(const IceModelVec2S &SeaLevel, double Offset) {
  m_floatation_level.copy_from(SeaLevel);
  this->fill2Level(SeaLevel, Offset);
}

void SeaLevelCC::fill2SeaLevel(const IceModelVec2S &SeaLevel, const IceModelVec2S &Offset) {
  m_floatation_level.copy_from(SeaLevel);
  this->fill2Level(SeaLevel, Offset);
}

void SeaLevelCC::labelMap_impl(unsigned int run_number,
                               std::vector<unsigned int> &i_vec,
                               std::vector<unsigned int> &j_vec,
                               std::vector<unsigned int> &parents,
                               std::vector<unsigned int> &lengths,
                               std::vector<bool> &isValidList) {
  (void) isValidList;
  //Set Sea-level to invalid value, where SL above bed, but run isolated from ocean
  IceModelVec::AccessList list{&m_floatation_level, &m_crop_mask};

  for(unsigned int k = 0; k <= run_number; ++k) {
    unsigned int label = trackParentRun(k, parents);
    for(unsigned int n = 0; n < lengths[k]; ++n) {
      if(label > 1) {
        const int i = i_vec[k] + n,
                  j = j_vec[k];
        m_floatation_level(i, j) = m_fill_value;
        m_crop_mask(i, j) = 1;
      }
    }
  }
}

void SeaLevelCC::prepare_mask_impl() {
  IceModelVec::AccessList list{&m_mask_run};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest   = (i == m_i_global_first),
         isEast   = (i == m_i_global_last),
         isSouth  = (j == m_j_global_first),
         isNorth  = (j == m_j_global_last),
         isMargin = (isWest or isEast or isSouth or isNorth);

    //Set sink at margin of computational domain
    m_mask_run(i, j) = isMargin ? 1 : 0;

  }
  m_mask_run.update_ghosts();
}

void SeaLevelCC::sea_level_2D(IceModelVec2S& result) const {
  this->get_floatation_level(result);
}

void SeaLevelCC::sl_crop_mask(IceModelVec2Int& result) const {
  result.copy_from(m_crop_mask);
}

}
