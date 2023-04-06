/* Copyright (C) 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "LakeAccumulator.hh"
#include "pism/util/error_handling.hh"

namespace pism {

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
        const int label = connected_components::trackParentRun(k, parents);
        for (unsigned int n = 0; n < len_vec[k]; ++n) {
          const int i = i_vec[k] + n;
          accumulator[label] += in_p0(i, j);
        }
      }

      //label result
      for (int k = 0; k <= m_run_number; ++k) {
        const int j = j_vec[k];
        const int label = connected_components::trackParentRun(k, parents);
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

} // end of namespace pism
