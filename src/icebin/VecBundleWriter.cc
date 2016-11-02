// See here for useful changes:
// https://github.com/pism/pism/commit/443050f30743d6c2ef431c53e87dc6eb19a73dfd

#include <base/util/PISMTime.hh>
#include <base/util/io/PIO.hh>
#include <base/util/io/io_helpers.hh>
#include <icebin/VecBundleWriter.hh>

using namespace pism;

namespace pism {
namespace icebin {

VecBundleWriter::VecBundleWriter(pism::IceGrid::Ptr _grid, std::string const &_fname,
                                 std::vector<pism::IceModelVec const *> &_vecs)
    : m_grid(_grid), fname(_fname), vecs(_vecs) {
}

void VecBundleWriter::init() {
  pism::PIO nc(m_grid->com, m_grid->ctx()->config()->get_string("output.format"),
               fname, PISM_READWRITE_MOVE);

  io::define_time(nc,
                  m_grid->ctx()->config()->get_string("time.dimension_name"),
                  m_grid->ctx()->time()->calendar(),
                  m_grid->ctx()->time()->CF_units_string(),
                  m_grid->ctx()->unit_system());

  for (pism::IceModelVec const *vec : vecs) {
    vec->define(nc, PISM_DOUBLE);
  }
}

/** Dump the value of the Vectors at curent PISM simulation time. */
void VecBundleWriter::write(double time_s) {
  pism::PIO nc(m_grid->com, m_grid->ctx()->config()->get_string("output.format"),
               fname.c_str(), PISM_READWRITE); // append to file

  io::append_time(nc, m_grid->ctx()->config()->get_string("time.dimension_name"), time_s);

  for (pism::IceModelVec const *vec : vecs) {
    vec->write(nc);
  }
}

} // end of namespace icebin
} // end of namespace pism
