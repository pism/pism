%{
  using pism::IceModelVec;
  using pism::IceModelVec2S;
  using pism::IceModelVec2Int;
  using pism::IceModelVec2V;
  using pism::IceModelVec3;
%}

/* disable methods that use regular pointers */
%ignore pism::Vars::add;
%ignore pism::Vars::is_available;
%ignore pism::Vars::get;
%ignore pism::Vars::get_2d_scalar;
%ignore pism::Vars::get_2d_vector;
%ignore pism::Vars::get_2d_mask;
%ignore pism::Vars::get_3d_scalar;
%ignore pism::Vars::keys;

/* replace with methods that use shared pointers */
%rename(add) pism::Vars::add_shared;
%rename(is_available) pism::Vars::is_available_shared;
%rename(get) pism::Vars::get_shared;
%rename(get_2d_scalar) pism::Vars::get_2d_scalar_shared;
%rename(get_2d_vector) pism::Vars::get_2d_vector_shared;
%rename(get_2d_mask) pism::Vars::get_2d_mask_shared;
%rename(get_3d_scalar) pism::Vars::get_3d_scalar_shared;
%rename(keys) pism::Vars::keys_shared;

/* remove Vars::remove */
%ignore pism::Vars::remove;

%extend pism::Vars
{
  %pythoncode
  {
    def __init__(self,*args):
      this = _cpp.new_Vars()
      try: self.this.append(this)
      except: self.this = this
      if len(args)>1:
        raise ValueError("Vars can only be constructed from nothing, an IceModelVec, or a list of such.")
      if len(args)==1:
        if isinstance(args[0],IceModelVec):
          self.add(args[0])
        else:
          # assume its a list of vecs
          for v in args[0]:
            self.add(v)
  }
}

%include "PISMVars.hh"
