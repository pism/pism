%{
/* Using directives needed to compile array::Array wrappers. */

#include "util/array/CellType.hh"
#include "util/array/Scalar.hh"
#include "util/array/Forcing.hh"
#include "util/array/Vector.hh"
#include "util/array/Array3D.hh"
#include "util/array/Staggered.hh"

using namespace pism;
%}

%shared_ptr(pism::PetscAccessible)
%shared_ptr(pism::array::Array)
%shared_ptr(pism::array::Scalar)
%shared_ptr(pism::array::Forcing)
%shared_ptr(pism::array::Vector)
%shared_ptr(pism::array::Vector1)
%shared_ptr(pism::array::Vector2)
%shared_ptr(pism::array::Staggered)
%shared_ptr(pism::array::Staggered1)
%shared_ptr(pism::array::Array3D)

%ignore pism::array::AccessScope::AccessScope(std::initializer_list<const PetscAccessible *>);

%ignore pism::array::Scalar::array;
%ignore pism::array::Vector::array;

%template(Range) std::array<double,2>;

%rename(_regrid) pism::array::Array::regrid;
%extend pism::array::Array
{
  %pythoncode "IceModelVec.py"
}

// Shenanigans to allow python indexing to get at array::Array entries.  I couldn't figure out a more
// elegant solution.
%extend pism::array::Scalar
{
    double getitem(int i, int j)
    {
        return (*($self))(i,j);
    }

    void setitem(int i, int j, double val)
    {
        (*($self))(i,j) = val;
    }

    %pythoncode "ArrayScalar.py"
};

%rename(__mult__) pism::Vector2d::operator*;
%rename(__add__) pism::Vector2d::operator+;
%ignore operator*(const double &a, const pism::Vector2d &v1);
%extend pism::Vector2d
{
  %pythoncode
  {
  def __lmul__(self,a):
    return self.__mul__(self,a)
  }
}

%extend pism::array::Vector
{
    Vector2d &getitem(int i, int j)
    {
        return (*($self))(i,j);
    }

    void setitem(int i, int j, Vector2d val)
    {
        (*($self))(i,j) = val;
    }

    void setitem(int i, int j, double u, double v)
    {
        (*($self))(i,j).u = u;
        (*($self))(i,j).v = v;
    }

    %pythoncode "IceModelVec2V.py"
};

%ignore pism::array::Array3D::operator();
%ignore pism::array::Array3D::get_column(int, int);
%ignore pism::array::Array3D::set_column(int, int, const double*);
%extend pism::array::Array3D
{

  double getitem(int i, int j, int k) const {
    return (*($self))(i,j,k);
  }

  void setitem(int i, int j, int k, double val) {
    (*($self))(i,j,k) = val;
  }

  std::vector<double> _get_column(int i, int j) const {
    size_t n = $self->levels().size();
    std::vector<double> result(n);
    const double *data = $self->get_column(i, j);
    for (size_t k = 0; k < n; ++k) {
      result[k] = data[k];
    }
    return result;
  }

  void set_column(int i, int j, const std::vector<double> &data) {
    assert(data.size() >= $self->levels().size());
    $self->set_column(i, j, data.data());
  }

    %pythoncode {
    def get_column(self, i, j):
          return self._get_column(i, j)

    def __getitem__(self,*args):
          i, j, k = args[0]
          return self.getitem(i, j, k)

    def __setitem__(self,*args):
        if(len(args)==2):
            i, j, k = args[0]
            self.setitem(i, j, k, args[1])
        else:
            raise ValueError("__setitem__ requires 2 arguments; received %d" % len(args));
    }
};

%ignore pism::array::Forcing::interp(int, int, double*);
%extend pism::array::Forcing
{
std::vector<double> interp(int i, int j) {
  std::vector<double> result;
  $self->interp(i, j, result);
  return result;
}
};

%ignore pism::stencils::Star::operator[];
%include "util/stencils.hh"
%template(DoubleStar) pism::stencils::Star<double>;

%include "util/array/Array.hh"

%include "util/array/Array2D.hh"
%shared_ptr(pism::array::Array2D<double>)
%ignore pism::array::Array2D< double >::array() const;
%template(_Array2Ddouble) pism::array::Array2D<double>;

%shared_ptr(pism::array::Scalar1)
%shared_ptr(pism::array::Scalar2)
%include "util/array/Scalar.hh"

%shared_ptr(pism::array::CellType)
%shared_ptr(pism::array::CellType1)
%shared_ptr(pism::array::CellType2)
%ignore pism::array::duplicate(const pism::array::CellType &);
%ignore pism::array::duplicate(const pism::array::CellType1 &);
%ignore pism::array::duplicate(const pism::array::CellType2 &);
%include "util/array/CellType.hh"

%include "util/array/Forcing.hh"

%shared_ptr(pism::array::Array2D<Vector2d>)
%ignore pism::array::Array2D< Vector2d >::array() const;
%template(_Array2DVector2) pism::array::Array2D<Vector2d>;
%include "util/array/Vector.hh"

%include "util/array/Array3D.hh"
%include "util/array/Staggered.hh"

%include "util/Vector2d.hh"
