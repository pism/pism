// See here for useful changes:
// https://github.com/pism/pism/commit/443050f30743d6c2ef431c53e87dc6eb19a73dfd

#include <pism/util/Time.hh>
#include <pism/util/io/File.hh>
#include <pism/util/io/io_helpers.hh>
#include <pism/icebin/VecBundleWriter.hh>
#include "pism/util/Context.hh"
#include "pism/util/array/Array.hh"
#include "pism/util/Grid.hh"

using namespace pism;

namespace pism {
namespace icebin {

VecBundleWriter::VecBundleWriter(std::shared_ptr<pism::Grid> _grid, std::string const &_fname,
                                 std::vector<pism::array::Array const *> &_vecs)
    : m_grid(_grid), fname(_fname), vecs(_vecs) {
}

void VecBundleWriter::init() {
  auto config = m_grid->ctx()->config();
  auto time = m_grid->ctx()->time();

  pism::File file(m_grid->com,
                  fname,
                  string_to_backend(config->get_string("output.format")),
                  io::PISM_READWRITE_MOVE,
                  m_grid->ctx()->pio_iosys_id());

  io::define_time(file,
                  config->get_string("time.dimension_name"),
                  time->calendar(),
                  time->units_string(),
                  m_grid->ctx()->unit_system());

  for (pism::array::Array const *vec : vecs) {
    vec->define(file, io::PISM_DOUBLE);
  }
}

/** Dump the value of the Vectors at curent PISM simulation time. */
void VecBundleWriter::write(double time_s) {
  auto config = m_grid->ctx()->config();

  pism::File file(m_grid->com,
                  fname,
                  string_to_backend(config->get_string("output.format")),
                  io::PISM_READWRITE, // append to file
                  m_grid->ctx()->pio_iosys_id());

  io::append_time(file, config->get_string("time.dimension_name"), time_s);

  for (pism::array::Array const *vec : vecs) {
    vec->write(file);
  }
}

} // end of namespace icebin
} // end of namespace pism
