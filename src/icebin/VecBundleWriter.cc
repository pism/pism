// See here for useful changes:
// https://github.com/pism/pism/commit/443050f30743d6c2ef431c53e87dc6eb19a73dfd

#include "pism/icebin/VecBundleWriter.hh"

#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Array.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Time.hh"

using namespace pism;

namespace pism {
namespace icebin {

VecBundleWriter::VecBundleWriter(std::shared_ptr<pism::Grid> _grid, std::string const &_fname,
                                 std::vector<pism::array::Array const *> &_vecs)
    : m_grid(_grid), fname(_fname), vecs(_vecs) {
}

void VecBundleWriter::init() {
  pism::File file(m_grid->com, fname,
                  string_to_backend(m_grid->ctx()->config()->get_string("output.format")),
                  io::PISM_READWRITE_MOVE);

  io::define_time(file, m_grid->ctx()->config()->get_string("time.dimension_name"),
                  m_grid->ctx()->time()->calendar(), m_grid->ctx()->time()->units_string(),
                  m_grid->ctx()->unit_system());

  for (const auto *vec : vecs) {
    vec->define(file, io::PISM_DOUBLE);
  }
}

/** Dump the value of the Vectors at curent PISM simulation time. */
void VecBundleWriter::write(double time_s) {
  pism::File nc(m_grid->com, fname,
                string_to_backend(m_grid->ctx()->config()->get_string("output.format")),
                io::PISM_READWRITE); // append to file

  io::append_time(nc, m_grid->ctx()->config()->get_string("time.dimension_name"), time_s);

  for (const auto *vec : vecs) {
    vec->write(nc);
  }
}

} // end of namespace icebin
} // end of namespace pism
