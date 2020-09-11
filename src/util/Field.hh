/* Copyright (C) 2014, 2015, 2019, 2020 PISM Authors
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

#ifndef _FIELD_H_
#define _FIELD_H_

namespace pism {

class Field {
protected:
  Field(IceGrid::ConstPtr g, const std::string &name);
public:
  virtual ~Field();
  // metadata
  std::string name() const;
  void set_name(const std::string &name, int component = 0);
  void set_units(const std::string &units, const std::string &external_units);

  SpatialVariableMetadata& metadata(unsigned int N = 0);
  const SpatialVariableMetadata& metadata(unsigned int N = 0) const;

  void set_attrs(const std::string &my_pism_intent, const std::string &my_long_name,
                 const std::string &my_units, const std::string &my_standard_name,
                 int component = 0);
  void read_attrs(const std::string &filename, int component = 0);

  // input and output
  void define(const File &nc, IO_Type output_datatype) const;

  void read(const std::string &filename, unsigned int time);
  void read(const File &nc, unsigned int time);

  void write(const std::string &filename, IO_Type nctype = PISM_DOUBLE) const;
  void write(const File &nc, IO_Type nctype = PISM_DOUBLE) const;

  void regrid(const std::string &filename, RegriddingFlag flag,
              double default_value = 0.0);
  void regrid(const File &nc, RegriddingFlag flag,
              double default_value = 0.0);

  // in-place modification
  void scale(double alpha);

  // low-level access
  Vec get_vec();
  PISMDM::Ptr get_dm() const;
  void copy_to_vec(PISMDM::Ptr destination_da, Vec destination) const;
  void copy_from_vec(Vec source);

  // point-wise access
  void begin_access() const;
  void end_access() const;

  // miscellaneous
  int modification_counter() const;
  int modified();
  void set_time_independent(bool value);

  // debugging
  void dump(const char filename[]) const;
protected:
  void get_dof(PISMDM::Ptr da_result, Vec result, unsigned int n,
               unsigned int count=1) const;
  void set_dof(PISMDM::Ptr da_source, Vec source, unsigned int n,
               unsigned int count=1);
  void check_array_indices(int i, int j, unsigned int k) const;

  std::string m_name;
  Vec m_v;
  IceGrid::ConstPtr m_grid;
  mutable void *m_array;
  //! distributed mesh manager (DM)
  PISMDM::Ptr m_da;
  unsigned int m_dof;

  //! stores metadata (NetCDF variable attributes)
  std::vector<SpatialVariableMetadata> m_metadata;
  std::string m_name;
private:
  mutable int m_access_counter;
  int m_modification_counter;

  // disable copy constructor and the assignment operator
  Field(const Field &);
  Field& operator=(const Field &);

public:
  //! Makes sure that we call begin_access() and end_access() for all accessed IceModelVecs.
  class AccessList {
  public:
    AccessList();
    AccessList(const IceModelVec &v);
    ~AccessList();
    void add(const IceModelVec &v);
  private:
    std::vector<const Field*> m_fields;
  };
};

class GhostedField : virtual public Field {
public:
  GhostedField(IceGrid::ConstPtr grid, const std::string &name, unsigned int stencil_width);
  virtual ~GhostedField();
  unsigned int stencil_width() const;
  void update_ghosts();
protected:
  unsigned int m_stencil_width;
  void scatter_to(Vec output) const;
};

class ScalarField : virtual public Field {
public:
  struct Range {
    double min, max;
  };
  ScalarField(IceGrid::ConstPtr grid, const std::string &name);
  virtual ~ScalarField();
  void set(double value);
  void shift(double amount);
  void squareroot();
  Range range() const;
};

class GhostedScalar2DField;
class Scalar2DField : virtual public ScalarField {
public:
  Scalar2DField(IceGrid::ConstPtr grid,  const std::string &name);
  virtual ~Scalar2DField();
  void copy_from(const Scalar2DField &input);
  void scatter_to_ghosted(GhostedScalar2DField &output) const;
  inline double& operator()(int i, int j);
  inline const double& operator()(int i, int j) const;
};

// finite differences
namespace FD {
inline double diff_x(const Scalar2DField &field, int i, int j) const;
inline double diff_y(const Scalar2DField &field, int i, int j) const;
inline double diff_x_p(const Scalar2DField &field, int i, int j) const;
inline double diff_y_p(const Scalar2DField &field, int i, int j) const;
}

class GhostedScalar3DField;
class Scalar3DField : virtual public ScalarField {
public:
  Scalar3DField(IceGrid::ConstPtr grid,  const std::string &name,
                const std::vector<double> &levels);
  virtual ~Scalar3DField();
  void copy_from(const Scalar3DField &input);
  void scatter_to_ghosted(GhostedScalar3DField &output) const;
  inline double& operator()(int i, int j, unsigned int k);
  inline const double& operator()(int i, int j, unsigned int k) const;
protected:
  std::vector<double> m_levels;
};

} // end of namespace pism

#endif /* _FIELD_H_ */
