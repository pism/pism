#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>
#include <pism/util/IceGrid.hh>
// --------------------------------

#include <string>
#include <vector>

namespace pism {

namespace array {
class Array;
} // end of namespace array

namespace icebin {


/** Sets up to easily write out a bundle of PISM variables to a file. */
class VecBundleWriter {
  pism::IceGrid::ConstPtr m_grid;
  std::string const fname;                     // Name of the file to write
  std::vector<pism::array::Array const *> vecs; // The vectors we will write

public:
  VecBundleWriter(std::shared_ptr<pism::IceGrid> grid, std::string const &_fname, std::vector<pism::array::Array const *> &_vecs);

  void init();

  /** Dump the value of the Vectors at curent PISM simulation time. */
  void write(double time_s);
};
} // end of namespace icebin
} // end of namespace pism
