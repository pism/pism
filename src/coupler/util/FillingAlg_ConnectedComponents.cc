#include <vector>

#include "pism/util/iceModelVec.hh"
#include "pism/util/pism_utilities.hh"

#include "FillingAlg_ConnectedComponents.hh"

namespace pism {

FillingAlgCC::FillingAlgCC(IceGrid::ConstPtr g, const double drho, const IceModelVec2S &bed, const IceModelVec2S &thk,
                           const bool check_validity, const double fill_value)
    : m_grid(g),
      m_drho(drho),
      m_i_local_first(m_grid->xs()),
      m_i_local_last(m_i_local_first + m_grid->xm() - 1),
      m_j_local_first(m_grid->ys()),
      m_j_local_last(m_j_local_first + m_grid->ym() - 1),
      m_i_global_first(0),
      m_i_global_last(m_grid->Mx() - 1),
      m_j_global_first(0),
      m_j_global_last(m_grid->My() - 1),
      m_bed(&bed),
      m_thk(&thk),
      m_check_validity(check_validity),
      m_fill_value(fill_value) {

  m_floatation_level.create(m_grid, "floatation_level", WITH_GHOSTS, 1);

  m_mask_run.create(m_grid, "mask_run", WITH_GHOSTS, 1);

  if (m_check_validity) {
    m_mask_validity.create(m_grid, "mask_validity", WITH_GHOSTS, 1);
  }
}

FillingAlgCC::~FillingAlgCC() {
  //empty
}

void FillingAlgCC::setLevel(double Level) {
  if (!m_level.was_created()) {
    m_level.create(m_grid, "current_level", WITHOUT_GHOSTS);
  }
  m_level.set(Level);
}

void FillingAlgCC::setLevel(const IceModelVec2S &Level) {
  if (!m_level.was_created()) {
    m_level.create(m_grid, "current_level", WITHOUT_GHOSTS);
  }
  m_level.copy_from(Level);
}

void FillingAlgCC::setOffset(double Offset) {
  if (!m_offset.was_created()) {
    m_offset.create(m_grid, "level_offset", WITHOUT_GHOSTS);
  }
  m_offset.set(Offset);
}

void FillingAlgCC::setOffset(const IceModelVec2S &Offset) {
  if (!m_offset.was_created()) {
    m_offset.create(m_grid, "level_offset", WITHOUT_GHOSTS);
  }
  m_offset.copy_from(Offset);
}

void FillingAlgCC::fill2Level(double Level) {
  this->setLevel(Level);
  this->fill2Level();
}

void FillingAlgCC::fill2Level(double Level, double Offset) {
  this->setLevel(Level);
  this->setOffset(Offset);
  this->fill2Level();
}

void FillingAlgCC::fill2Level(const IceModelVec2S &Level) {
  this->setLevel(Level);
  this->fill2Level();
}

void FillingAlgCC::fill2Level(const IceModelVec2S &Level, double Offset) {
  this->setLevel(Level);
  this->setOffset(Offset);
  this->fill2Level();
}

void FillingAlgCC::fill2Level(const IceModelVec2S &Level, const IceModelVec2S &Offset) {
  this->setLevel(Level);
  this->setOffset(Offset);
  this->fill2Level();
}

void FillingAlgCC::fill2Level() {

  if (!m_level.was_created()) {
    this->setLevel(m_fill_value);
  }

  if (!m_offset.was_created()) {
    this->setOffset(m_fill_value);
  }

  unsigned int max_items = 2 * m_grid->ym();

  std::vector<unsigned int> parents(max_items), lengths(max_items), j_vec(max_items), i_vec(max_items);
  std::vector<bool> isValidList(max_items);

  //Initialize vectors
  for (unsigned int k = 0; k < 2; ++k) {
    parents[k]     = 0;
    lengths[k]     = 0;
    j_vec[k]       = 0;
    i_vec[k]       = 0;
    isValidList[k] = true;
  }

  //run_number = 0 reserved for 'Background' Pixels,
  //run_number = 1 reserved for overflowing Pixel (or Sink)
  unsigned int run_number = 1;
  IceModelVec::AccessList list{ m_bed, m_thk, &m_mask_run, &m_level, &m_offset };
  if (m_check_validity) {
    list.add(m_mask_validity);
  }

  //Assign Pixels to runs
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double bed = (*m_bed)(i, j), thk = (*m_thk)(i, j),
           level = m_level(i, j), offset = m_offset(i, j);
    StarStencil<int> mask_star = m_mask_run.int_star(i, j);

    if (ForegroundCond(bed, thk, mask_star.ij, level, offset)) {

      bool isValid = true;
      if (m_check_validity) {
        isValid = (m_mask_validity(i, j) > 0);
      }

      checkForegroundPixel((unsigned int)i, (unsigned int)j, SinkCond(mask_star.ij), isValid, mask_star,
                           run_number, i_vec, j_vec, parents, lengths, isValidList);

      m_mask_run(i, j) = run_number;

      //resize vectors if 'max_items' are exceeded
      if ((run_number + 1) >= max_items) {
        max_items += m_grid->ym();
        parents.resize(max_items);
        lengths.resize(max_items);
        j_vec.resize(max_items);
        i_vec.resize(max_items);
        isValidList.resize(max_items);
      }
    }
  }

  labelMask(run_number, i_vec, j_vec, parents, lengths, isValidList);

  bool updateAtBoundaries = true;

  //iteratively adapt fields amongst processor domains
  while (updateAtBoundaries) {
    updateAtBoundaries = updateRunsAtBoundaries(parents, isValidList);

    labelMask(run_number, i_vec, j_vec, parents, lengths, isValidList);
  }

  labelMap(run_number, i_vec, j_vec, parents, lengths, isValidList);
}

void FillingAlgCC::checkForegroundPixel(unsigned int i, unsigned int j, bool isSink, bool isValid,
                                        StarStencil<int> mask_star, unsigned int &run_number,
                                        std::vector<unsigned int> &i_vec, std::vector<unsigned int> &j_vec,
                                        std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths,
                                        std::vector<bool> &isValidList) {
  bool isWest = (i <= m_i_local_first), isSouth = (j <= m_j_local_first);

  if (not isWest and (mask_star.w > 0)) {
    // west neighbor is also foreground: continue the run
    lengths[run_number] += 1;
    if (isSink) {
      setRunSink(run_number, parents);
    }
    if (m_check_validity and isValid) {
      setRunValid(run_number, parents, isValidList);
    }
  } else {
    //west neighbor is a background pixel (or this is westmost column): start a new run
    unsigned int parent;
    if (isSink) {
      parent = 1;
    } else if (not isSouth and (mask_star.s > 0)) {
      //check the pixel south and set the parent
      parent = (unsigned int)mask_star.s;
    } else {
      parent = 0;
    }

    run_number += 1;

    i_vec[run_number]       = i;
    j_vec[run_number]       = j;
    parents[run_number]     = parent;
    lengths[run_number]     = 1;
    isValidList[run_number] = isValid;
  }

  if (not isSouth and (mask_star.s > 0)) {
    run_union(parents, (unsigned int)mask_star.s, run_number);
    if (isValid and m_check_validity) {
      setRunValid(run_number, parents, isValidList);
    }
  }
}

bool FillingAlgCC::updateRunsAtBoundaries(std::vector<unsigned int> &parents, std::vector<bool> &isValidList) {
  bool changed = false;
  IceModelVec::AccessList list{ &m_mask_run };
  m_mask_run.update_ghosts();

  if (m_check_validity) {
    list.add(m_mask_validity);
    m_mask_validity.update_ghosts();
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest   = ((i == m_i_local_first) and not(i == m_i_global_first)),
         isEast   = ((i == m_i_local_last) and not(i == m_i_global_last)),
         isSouth  = ((j == m_j_local_first) and not(j == m_j_global_first)),
         isNorth  = ((j == m_j_local_last) and not(j == m_j_global_last)),
         isMargin = (isWest or isEast or isSouth or isNorth);

    if (isMargin) {
      unsigned int run = m_mask_run.as_int(i, j);
      if (run > 1) {
        //Lake at inner boundary
        StarStencil<int> mask_star = m_mask_run.int_star(i, j);
        bool WestSink = (isWest and (mask_star.w == 1)), EastSink = (isEast and (mask_star.e == 1)),
             SouthSink = (isSouth and (mask_star.s == 1)), NorthSink = (isNorth and (mask_star.n == 1));

        if (WestSink or EastSink or SouthSink or NorthSink) {
          //Lake on other side overflowing
          parents[run] = 1;
          changed      = true;
        } else if (m_check_validity and not isValidList[run]) {
          //Lake on other side not overflowing and
          //Lake at this side is not labeled as invalid
          StarStencil<int> mask_isValid_star = m_mask_validity.int_star(i, j);

          bool WestValid  = (isWest  and (mask_isValid_star.w == 1)),
               EastValid  = (isEast  and (mask_isValid_star.e == 1)),
               SouthValid = (isSouth and (mask_isValid_star.s == 1)),
               NorthValid = (isNorth and (mask_isValid_star.n == 1));

          if (WestValid or EastValid or SouthValid or NorthValid) {
            //Lake at other side is not completely covered with ice
            isValidList[run] = true;
            changed          = true;
          }
        }
      }
    }
  }

  return (GlobalOr(m_grid->com, changed));
}

void FillingAlgCC::labelMask_impl(unsigned int run_number, std::vector<unsigned int> &i_vec,
                                  std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents,
                                  std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) {
  int label_valid = 1;

  IceModelVec::AccessList list{ &m_mask_run };
  if (m_check_validity) {
    list.add(m_mask_validity);
  }
  for (unsigned int k = 0; k <= run_number; ++k) {
    unsigned int label = trackParentRun(k, parents);
    if (m_check_validity) {
      label_valid = isValidList[label] ? 1 : 0;
    }
    for (unsigned int n = 0; n < lengths[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      m_mask_run(i, j) = label;
      if (m_check_validity) {
        m_mask_validity(i, j) = label_valid;
      }
    }
  }
}

void FillingAlgCC::run_union(std::vector<unsigned int> &parents, unsigned int run1, unsigned int run2) {
  if ((parents[run1] == run2) or (parents[run2] == run1)) {
    return;
  }

  run1 = trackParentRun(run1, parents);

  run2 = trackParentRun(run2, parents);

  if (run1 > run2) {
    parents[run1] = run2;
  } else if (run1 < run2) {
    parents[run2] = run1;
  }
}

unsigned int FillingAlgCC::trackParentRun(unsigned int run, std::vector<unsigned int> &parents) {
  while (parents[run] != 0) {
    run = parents[run];
  }
  return run;
}

void FillingAlgCC::setRunSink(unsigned int run, std::vector<unsigned int> &parents) {
  if ((run == 0) or (run == 1)) {
    return;
  }

  run = trackParentRun(run, parents);
  if (run != 1) {
    parents[run] = 1;
  }
}

void FillingAlgCC::setRunValid(unsigned int run, std::vector<unsigned int> &parents, std::vector<bool> &isValidList) {
  if ((run == 0) or (run == 1)) {
    return;
  }

  run = trackParentRun(run, parents);
  if (run != 1) {
    isValidList[run] = true;
  }
}

bool FillingAlgCC::SinkCond_impl(int mask) const {
  return (mask == 1);
}

bool FillingAlgCC::ForegroundCond_impl(double bed, double thk, int mask, double Level, double Offset) const {
  if (mask > 0) {
    return true;
  }

  if (Level == m_fill_value) {
    return true;
  }

  double level = Level;
  if (Offset != m_fill_value) {
    level += Offset;
  }

  return ((bed + (m_drho * thk)) < level);
}

void FillingAlgCC::get_floatation_level(IceModelVec2S &result) const {
  result.copy_from(m_floatation_level);
}

void FillingAlgCC::prepare_mask() {
  this->prepare_mask_impl();
}

bool FillingAlgCC::SinkCond(int mask) const {
  return this->SinkCond_impl(mask);
}

bool FillingAlgCC::ForegroundCond(double bed, double thk, int mask, double Level, double Offset) const {
  return this->ForegroundCond_impl(bed, thk, mask, Level, Offset);
}

void FillingAlgCC::labelMask(unsigned int run_number, std::vector<unsigned int> &i_vec,
                             std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents,
                             std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) {
  this->labelMask_impl(run_number, i_vec, j_vec, parents, lengths, isValidList);
}

void FillingAlgCC::labelMap(unsigned int run_number, std::vector<unsigned int> &i_vec,
                            std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents,
                            std::vector<unsigned int> &lengths, std::vector<bool> &isValidList) {
  this->labelMap_impl(run_number, i_vec, j_vec, parents, lengths, isValidList);
}

} //namespace pism
