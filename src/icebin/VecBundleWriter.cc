// See here for useful changes:
// https://github.com/pism/pism/commit/443050f30743d6c2ef431c53e87dc6eb19a73dfd

#include "pism/icebin/VecBundleWriter.hh"

#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Array.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"

using namespace pism;

namespace pism {
namespace icebin {

VecBundleWriter::VecBundleWriter(std::shared_ptr<pism::Grid> _grid, std::string const &_fname,
                                 std::vector<pism::array::Array const *> &_vecs)
    : m_grid(_grid),
      fname(_fname),
      vecs(_vecs),
      output_writer(new SynchronousOutputWriter(m_grid->com, *m_grid->ctx()->config())) {
}

void VecBundleWriter::init() {
  pism::OutputFile file(output_writer, fname);

  file.define_variable(m_grid->ctx()->time()->metadata());
  
  for (const auto *vec : vecs) {
    for (const auto &var : vec->all_metadata()) {
      file.define_variable(var);
    }
  }
}

/** Dump the value of the Vectors at curent PISM simulation time. */
void VecBundleWriter::write(double time_s) {
  pism::OutputFile file(output_writer, fname);

  file.append_time(time_s);
  for (const auto *vec : vecs) {
    vec->write(file);
  }
}

} // end of namespace icebin
} // end of namespace pism
