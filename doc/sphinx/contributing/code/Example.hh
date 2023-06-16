#ifndef _PO_EXAMPLE_H_
#define _PO_EXAMPLE_H_

#include "coupler/PISMOcean.hh"
#include "base/util/iceModelVec2T.hh"
#include <memory>

namespace pism {
namespace ocean {
//! \brief An example ocean model illustrating the use of `array::Forcing`.
class Example : public OceanModel {
public:
  Example(std::shared_ptr<const Grid> grid);
  virtual ~Example();
protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(double my_t, double my_dt);
  virtual void init_impl();
  virtual void sea_level_elevation_impl(double &result) const;
  virtual void shelf_base_temperature_impl(array::Scalar &result) const;
  virtual void shelf_base_mass_flux_impl(array::Scalar &result) const;
protected:
  array::Forcing m_shelf_melt_rate;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _PO_EXAMPLE_H_ */
