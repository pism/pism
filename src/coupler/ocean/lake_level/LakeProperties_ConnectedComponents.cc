#include <vector>

#include "pism/util/iceModelVec.hh"
#include "pism/util/pism_utilities.hh"

#include "LakeProperties_ConnectedComponents.hh"

namespace pism {

LakePropertiesCC::LakePropertiesCC(IceGrid::ConstPtr g, const double drho, const double fill_value, const IceModelVec2S *target_level, const IceModelVec2S *lake_level, const IceModelVec2S *floating_thresh)
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
      m_target_level(target_level),
      m_current_level(lake_level),
      m_floating_threshold_level(floating_thresh),
      m_fill_value(fill_value) {

  m_mask_run.create(m_grid, "mask_run", WITH_GHOSTS, 1);
  m_mask_run.set(0);
  m_min_ll_tmp.create(m_grid, "min_ll_tmp", WITH_GHOSTS, 1);
  m_min_ll_tmp.set(m_fill_value);
  m_max_ll_tmp.create(m_grid, "max_ll_tmp", WITH_GHOSTS, 1);
  m_max_ll_tmp.set(m_fill_value);
  m_min_float_tmp.create(m_grid, "min_float_tmp", WITH_GHOSTS, 1);
  m_min_float_tmp.set(m_fill_value);
}

LakePropertiesCC::~LakePropertiesCC() {
  //empty
}

void LakePropertiesCC::getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level, IceModelVec2S &min_float_level) {
  unsigned int max_items = 2 * m_grid->ym();

  std::vector<unsigned int> parents(max_items), lengths(max_items), j_vec(max_items), i_vec(max_items);
  std::vector<double> minlevel_list(max_items), maxlevel_list(max_items), minfloatlevel_list(max_items);

  //Initialize vectors
  parents[0]   = 0;
  lengths[0]   = 0;
  j_vec[0]     = 0;
  i_vec[0]     = 0;
  minlevel_list[0]      = m_fill_value;
  maxlevel_list[0]      = m_fill_value;
  minfloatlevel_list[0] = m_fill_value;

  //run_number = 0 reserved for 'Background' Pixels,
  unsigned int run_number = 0;
  IceModelVec::AccessList list{ m_target_level, m_current_level, m_floating_threshold_level, &m_mask_run };

  //Assign Pixels to runs
  for (Points p(*m_grid); p; p.next()) {
    const unsigned int i = p.i(), j = p.j();

    const double target_ij       = (*m_target_level)(i, j),
                 current_ij      = (*m_current_level)(i, j),
                 float_thresh_ij = (*m_floating_threshold_level)(i, j);
    StarStencil<int> mask_star = m_mask_run.int_star(i, j);

    if (ForegroundCond(target_ij, current_ij)) {

      checkForegroundPixel(i, j, mask_star, run_number, current_ij, float_thresh_ij, i_vec, j_vec, parents, lengths, minlevel_list, maxlevel_list, minfloatlevel_list);

      m_mask_run(i, j) = run_number;

      //resize vectors if 'max_items' are exceeded
      if ((run_number + 1) >= max_items) {
        max_items += m_grid->ym();
        parents.resize(max_items);
        lengths.resize(max_items);
        j_vec.resize(max_items);
        i_vec.resize(max_items);
        minlevel_list.resize(max_items);
        maxlevel_list.resize(max_items);
        minfloatlevel_list.resize(max_items);
      }
    }
  }

  labelMask(run_number, i_vec, j_vec, parents, lengths, minlevel_list, maxlevel_list, minfloatlevel_list);

  bool updateAtBoundaries = true;

  //iteratively adapt fields amongst processor domains
  while (updateAtBoundaries) {
    updateAtBoundaries = updateRunsAtBoundaries(parents, minlevel_list, maxlevel_list, minfloatlevel_list);

    labelMask(run_number, i_vec, j_vec, parents, lengths, minlevel_list, maxlevel_list, minfloatlevel_list);
  }

//   labelMap(Level, run_number, i_vec, j_vec, parents, lengths, minlevel_list, maxlevel_list, minfloatlevel_list);

min_level.copy_from(m_min_ll_tmp);
max_level.copy_from(m_max_ll_tmp);
min_float_level.copy_from(m_min_float_tmp);
}

void LakePropertiesCC::checkForegroundPixel(unsigned int i, unsigned int j, StarStencil<int> mask_star, unsigned int &run_number, const double lake_level, const double float_thresh, std::vector<unsigned int> &i_vec, std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list) {
  bool isWest = (i <= m_i_local_first), isSouth = (j <= m_j_local_first);

  if (not isWest and (mask_star.w > 0)) {
    // west neighbor is also foreground: continue the run
    lengths[run_number] += 1;
    setRunMinLevel(run_number, parents, minlevel_list, lake_level);
    setRunMaxLevel(run_number, parents, maxlevel_list, lake_level);
    setRunMinLevel(run_number, parents, minfloatlevel_list, float_thresh);
  } else {
    //west neighbor is a background pixel (or this is westmost column): start a new run
    unsigned int parent;
    if (not isSouth and (mask_star.s > 0)) {
      //check the pixel south and set the parent
      parent = (unsigned int)mask_star.s;
    } else {
      parent = 0;
    }

    run_number += 1;

    i_vec[run_number]     = i;
    j_vec[run_number]     = j;
    parents[run_number]   = parent;
    lengths[run_number]   = 1;
    minlevel_list[run_number]      = lake_level;
    maxlevel_list[run_number]      = lake_level;
    minfloatlevel_list[run_number] = float_thresh;
  }

  if (not isSouth and (mask_star.s > 0)) {
    run_union(parents, (unsigned int)mask_star.s, run_number);
  }
}

bool LakePropertiesCC::updateRunsAtBoundaries(std::vector<unsigned int> &parents, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list) {
  bool changed = false;
  IceModelVec::AccessList list{ &m_mask_run, &m_min_ll_tmp, &m_max_ll_tmp, &m_min_float_tmp };
  m_mask_run.update_ghosts();
  m_min_ll_tmp.update_ghosts();
  m_max_ll_tmp.update_ghosts();
  m_min_float_tmp.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    bool isWest   = ((i == m_i_local_first) and not(i == m_i_global_first)),
         isEast   = ((i == m_i_local_last) and not(i == m_i_global_last)),
         isSouth  = ((j == m_j_local_first) and not(j == m_j_global_first)),
         isNorth  = ((j == m_j_local_last) and not(j == m_j_global_last)),
         isMargin = (isWest or isEast or isSouth or isNorth);

    if (isMargin) {
      unsigned int run = m_mask_run.as_int(i, j);
      if (run > 0) {
        StarStencil<double> min_star      = m_min_ll_tmp.star(i, j),
                            max_star      = m_max_ll_tmp.star(i, j),
                            minfloat_star = m_min_float_tmp.star(i, j);

        double min_level = min_star.ij,
              max_level = max_star.ij,
              min_float = minfloat_star.ij;

        if (isWest) {
          if (((min_star.w != m_fill_value) and (min_star.w < min_level)) or (min_level == m_fill_value)) {
            min_level = min_star.w;
          }
          if (((max_star.w != m_fill_value) and (max_star.w > max_level)) or (max_level == m_fill_value)) {
            max_level = max_star.w;
          }
          if (((minfloat_star.w != m_fill_value) and (minfloat_star.w < min_float)) or (min_float == m_fill_value)) {
            min_float = minfloat_star.w;
          }
        }
        if (isNorth) {
          if (((min_star.n != m_fill_value) and (min_star.n < min_level)) or (min_level == m_fill_value)) {
            min_level = min_star.n;
          }
          if (((max_star.n != m_fill_value) and (max_star.n > max_level)) or (max_level == m_fill_value)) {
            max_level = max_star.n;
          }
          if (((minfloat_star.n != m_fill_value) and (minfloat_star.n < min_float)) or (min_float == m_fill_value)) {
            min_float = minfloat_star.n;
          }
        }
        if (isEast) {
          if (((min_star.e != m_fill_value) and (min_star.e < min_level)) or (min_level == m_fill_value)) {
            min_level = min_star.e;
          }
          if (((max_star.e != m_fill_value) and (max_star.e > max_level)) or (max_level == m_fill_value)) {
            max_level = max_star.e;
          }
          if (((minfloat_star.e != m_fill_value) and (minfloat_star.e < min_float)) or (min_float == m_fill_value)) {
            min_float = minfloat_star.e;
          }
        }
        if (isSouth) {
          if (((min_star.s != m_fill_value) and (min_star.s < min_level)) or (min_level == m_fill_value)) {
            min_level = min_star.s;
          }
          if (((max_star.s != m_fill_value) and (max_star.s > max_level)) or (max_level == m_fill_value)) {
            max_level = max_star.s;
          }
          if (((minfloat_star.s != m_fill_value) and (minfloat_star.s < min_float)) or (min_float == m_fill_value)) {
            min_float = minfloat_star.s;
          }
        }
        if (min_level != min_star.ij) {
          setRunMinLevel(run, parents, minlevel_list, min_level);
          changed = true;
        }
        if (max_level != max_star.ij) {
          setRunMaxLevel(run, parents, maxlevel_list, max_level);
          changed = true;
        }
        if (min_float != minfloat_star.ij) {
          setRunMinLevel(run, parents, minfloatlevel_list, min_float);
          changed = true;
        }
      }
    }
  }
  return (GlobalOr(m_grid->com, changed));
}

void LakePropertiesCC::labelMask_impl(unsigned int run_number, std::vector<unsigned int> &i_vec, std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list) {
  IceModelVec::AccessList list{ &m_mask_run, &m_min_ll_tmp, &m_max_ll_tmp, &m_min_float_tmp };

  for (unsigned int k = 0; k <= run_number; ++k) {
    unsigned int label = trackParentRun(k, parents);
    const double min_ll_label        = minlevel_list[label],
                 max_ll_label        = maxlevel_list[label],
                 minfloatlevel_label = minfloatlevel_list[label];
    for (unsigned int n = 0; n < lengths[k]; ++n) {
      const int i = i_vec[k] + n, j = j_vec[k];
      m_mask_run(i, j)      = label;
      m_min_ll_tmp(i, j)    = min_ll_label;
      m_max_ll_tmp(i, j)    = max_ll_label;
      m_min_float_tmp(i, j) = minfloatlevel_label;
    }
  }
}

void LakePropertiesCC::run_union(std::vector<unsigned int> &parents, unsigned int run1, unsigned int run2) {
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

unsigned int LakePropertiesCC::trackParentRun(unsigned int run, std::vector<unsigned int> &parents) {
  while (parents[run] != 0) {
    run = parents[run];
  }
  return run;
}

void LakePropertiesCC::setRunMinLevel(unsigned int run, std::vector<unsigned int> &parents, std::vector<double> &run_level_list, double level) {
  if (run == 0) {
    return;
  }

  run = trackParentRun(run, parents);
  if (level != m_fill_value) {
    if (run_level_list[run] != m_fill_value) {
      run_level_list[run] = std::min(level, run_level_list[run]);
    } else {
      run_level_list[run] = level;
    }
  }
}

void LakePropertiesCC::setRunMaxLevel(unsigned int run, std::vector<unsigned int> &parents, std::vector<double> &run_level_list, double level) {
  if (run == 0) {
    return;
  }

  run = trackParentRun(run, parents);
  if (level != m_fill_value) {
    if (run_level_list[run] != m_fill_value) {
      run_level_list[run] = std::max(level, run_level_list[run]);
    } else {
      run_level_list[run] = level;
    }
  }
}

bool LakePropertiesCC::ForegroundCond_impl(double target_ll, double current_ll) const {
  return ((target_ll != m_fill_value) or (current_ll != m_fill_value));
}

bool LakePropertiesCC::ForegroundCond(double target_ll, double current_ll) const {
  return this->ForegroundCond_impl(target_ll, current_ll);
}

void LakePropertiesCC::labelMask(unsigned int run_number, std::vector<unsigned int> &rows, std::vector<unsigned int> &columns, std::vector<unsigned int> &parents, std::vector<unsigned int> &lengths, std::vector<double> &minlevel_list, std::vector<double> &maxlevel_list, std::vector<double> &minfloatlevel_list) {
  this->labelMask_impl(run_number, rows, columns, parents, lengths, minlevel_list, maxlevel_list, minfloatlevel_list);
}

// void LakePropertiesCC::labelMap(double Level, unsigned int run_number, std::vector<unsigned int> &i_vec,
//                             std::vector<unsigned int> &j_vec, std::vector<unsigned int> &parents,
//                             std::vector<unsigned int> &lengths, std::vector<bool> &isIceFree) {
//   this->labelMap_impl(Level, run_number, i_vec, j_vec, parents, lengths, isIceFree);
// }

} //namespace pism
