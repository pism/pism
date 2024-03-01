#ifndef PISM_YACINTERPOLATION_H
#define PISM_YACINTERPOLATION_H

#include <memory>
#include <string>

#include "pism/util/io/IO_Flags.hh"
#include "pism/util/Units.hh"
#include "pism/util/projection.hh"

namespace pism {
class Grid;
class File;

namespace array {
class Scalar;
}
} // namespace pism

class YACInterpolation {
public:
  YACInterpolation(const pism::Grid &grid,
                   const pism::File &file,
                   const std::string &variable_name);
  ~YACInterpolation();

  void regrid(const pism::File &file, pism::io::Default default_value,
              pism::array::Scalar &target) const;

  static std::string grid_name(const pism::File &file, const std::string &variable_name,
                      pism::units::System::Ptr sys);
private:
  double interpolate(const pism::array::Scalar &source,
                     pism::array::Scalar &target) const;

  static int interpolation_coarse_to_fine(double missing_value);
  static int interpolation_fine_to_coarse(double missing_value);

  static int define_field(int component_id, const pism::Grid &pism_grid,
                          const std::string &name);
  static int define_grid(const pism::Grid &grid, const std::string &grid_name,
                         const std::string &projection);

  int m_instance_id;
  int m_source_field_id;
  int m_target_field_id;

  std::shared_ptr<pism::array::Scalar> m_buffer;
};

pism::MappingInfo mapping(const pism::File &file,
                          pism::units::System::Ptr sys);

#endif /* PISM_YACINTERPOLATION_H */
