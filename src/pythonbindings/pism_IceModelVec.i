%{
/* Using directives needed to compile IceModelVec wrappers. */
#include "util/IceModelVec2CellType.hh"
#include "util/iceModelVec2T.hh"

using namespace pism;
%}

%shared_ptr(pism::PetscAccessible)
%shared_ptr(pism::IceModelVec)
%shared_ptr(pism::IceModelVec2S)
%shared_ptr(pism::IceModelVec2T)
%shared_ptr(pism::IceModelVec2V)
%shared_ptr(pism::IceModelVec2Int)
%shared_ptr(pism::IceModelVec2CellType)
%shared_ptr(pism::IceModelVec2Stag)
%shared_ptr(pism::IceModelVec3)

%ignore pism::AccessList::AccessList(std::initializer_list<const PetscAccessible *>);

%ignore pism::IceModelVec2S::array;
%ignore pism::IceModelVec2V::array;

%template(Range) std::array<double,2>;

%rename(_regrid) pism::IceModelVec::regrid;
%extend pism::IceModelVec
{
  %pythoncode "IceModelVec.py"
}

// Shenanigans to allow python indexing to get at IceModelVec entries.  I couldn't figure out a more
// elegant solution.
%extend pism::IceModelVec2S
{
    double getitem(int i, int j)
    {
        return (*($self))(i,j);
    }

    void setitem(int i, int j, double val)
    {
        (*($self))(i,j) = val;
    }

    %pythoncode "IceModelVec2S.py"
};

%rename(__mult__) pism::Vector2::operator*;
%rename(__add__) pism::Vector2::operator+;
%ignore operator*(const double &a, const pism::Vector2 &v1);
%extend pism::Vector2
{
  %pythoncode
  {
  def __lmul__(self,a):
    return self.__mul__(self,a)
  }
}

%extend pism::IceModelVec2V
{
    Vector2 &getitem(int i, int j)
    {
        return (*($self))(i,j);
    }

    void setitem(int i, int j, Vector2 val)
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

%ignore pism::IceModelVec3::operator();
%ignore pism::IceModelVec3::get_column(int, int);
%ignore pism::IceModelVec3::set_column(int, int, const double*);
%extend pism::IceModelVec3
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

%ignore pism::IceModelVec2T::interp(int, int, double*);
%extend pism::IceModelVec2T
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

%include "util/iceModelVec.hh"
%include "util/IceModelVec2.hh"

%shared_ptr(pism::IceModelVec2<Vector2>)
%ignore pism::IceModelVec2< Vector2 >::array() const;
%template(_IceModelVec2Vector2) pism::IceModelVec2<Vector2>;

%include "util/IceModelVec2V.hh"
%include "util/IceModelVec2CellType.hh"
%include "util/iceModelVec2T.hh"
%include "util/Vector2.hh"
