%{
/* Using directives needed to compile IceModelVec wrappers. */
#include "util/IceModelVec2CellType.hh"
#include "util/iceModelVec2T.hh"

using namespace pism;
%}

%shared_ptr(pism::PetscAccessible)
%shared_ptr(pism::IceModelVec)
%shared_ptr(pism::IceModelVec2)
%shared_ptr(pism::IceModelVec2S)
%shared_ptr(pism::IceModelVec2T)
%shared_ptr(pism::IceModelVec2V)
%shared_ptr(pism::IceModelVec2Int)
%shared_ptr(pism::IceModelVec2CellType)
%shared_ptr(pism::IceModelVec2Stag)
%shared_ptr(pism::IceModelVec3D)
%shared_ptr(pism::IceModelVec3)

%ignore pism::AccessList::AccessList(std::initializer_list<const PetscAccessible *>);

%ignore pism::IceModelVec2S::get_array;
%ignore pism::IceModelVec2V::get_array;

%rename(_regrid) pism::IceModelVec::regrid;
%extend pism::IceModelVec
{
  %pythoncode "IceModelVec.py"
}

// We also make the same fix for IceModelVec2's.
%rename(_regrid) pism::IceModelVec2::regrid;
%extend pism::IceModelVec2
{
  %pythoncode {
    def regrid(self,filename,critical=False,default_value=0.0):
      if critical == True:
        flag = CRITICAL
      else:
        flag = OPTIONAL
      self._regrid(filename, flag, default_value)
  }
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

%ignore pism::IceModelVec3D::getInternalColumn(int,int) const;
%ignore pism::IceModelVec3D::operator();
%extend pism::IceModelVec3D
{

  double getitem(int i, int j, int k)
  {
      return (*($self))(i,j,k);
  }

  void setitem(int i, int j, int k, double val)
  {
      (*($self))(i,j,k) = val;
  }


    %pythoncode {
    def __getitem__(self,*args):
        return self.getitem(args[0][0],args[0][1],args[0][2])

    def __setitem__(self,*args):
        if(len(args)==2):
            self.setitem(args[0][0],args[0][1],args[0][2],args[1])
        else:
            raise ValueError("__setitem__ requires 2 arguments; received %d" % len(args));
    }
};

%ignore pism::StarStencil::operator[];
%include "util/iceModelVec.hh"
%include "util/IceModelVec2CellType.hh"
%include "util/iceModelVec2T.hh"
%include "util/Vector2.hh"
