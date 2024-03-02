#ifndef PISM_YACINTERPOLATION_H
#define PISM_YACINTERPOLATION_H

#include <memory>
#include <string>

#include "pism/util/io/IO_Flags.hh"
#include "pism/util/Units.hh"

namespace pism {
class Grid;
class File;

namespace array {
class Scalar;
}

class YACInterpolation {
public:
  YACInterpolation(const Grid &target_grid, const File &input_file,
                   const std::string &variable_name);
  ~YACInterpolation();

  void regrid(const File &file, io::Default default_value, array::Scalar &target) const;

  static std::string grid_name(const File &file, const std::string &variable_name,
                               units::System::Ptr sys);

private:
  double interpolate(const array::Scalar &source, array::Scalar &target) const;

  static int interpolation_coarse_to_fine(double missing_value);
  static int interpolation_fine_to_coarse(double missing_value);

  static int define_field(int component_id, const Grid &pism_grid, const std::string &name);
  static int define_grid(const Grid &grid, const std::string &grid_name,
                         const std::string &projection);

  int m_instance_id;
  int m_source_field_id;
  int m_target_field_id;

  std::shared_ptr<array::Scalar> m_buffer;
};

} // namespace pism

#endif /* PISM_YACINTERPOLATION_H */
