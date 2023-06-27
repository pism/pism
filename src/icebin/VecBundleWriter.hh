#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <base/iceModel.hh>
#include <base/util/IceGrid.hh>
#include <base/util/io/PISMNCFile.hh>
#include <petsc.h>
// --------------------------------

#include <string>
#include <vector>

namespace pism {
namespace icebin {


/** Sets up to easily write out a bundle of PISM variables to a file. */
class VecBundleWriter {
  pism::IceGrid::ConstPtr m_grid;
  std::string const fname;                     // Name of the file to write
  std::vector<pism::IceModelVec const *> vecs; // The vectors we will write

public:
  VecBundleWriter(pism::IceGrid::Ptr grid, std::string const &_fname, std::vector<pism::IceModelVec const *> &_vecs);

  void init();

  /** Dump the value of the Vectors at curent PISM simulation time. */
  void write(double time_s);
};
} // end of namespace icebin
} // end of namespace pism
