%{
  #include "util/Vars.hh"
%}

/* disable methods that use regular pointers */
%ignore pism::Vars::add;
%ignore pism::Vars::is_available;
%ignore pism::Vars::get;
%ignore pism::Vars::get_2d_scalar;
%ignore pism::Vars::get_2d_vector;
%ignore pism::Vars::get_2d_mask;
%ignore pism::Vars::get_2d_cell_type;
%ignore pism::Vars::get_3d_scalar;
%ignore pism::Vars::keys;

/* replace with methods that use shared pointers */
%rename(add) pism::Vars::add_shared;
%rename(is_available) pism::Vars::is_available_shared;
%rename(_get) pism::Vars::get_shared;
%rename(get_2d_scalar) pism::Vars::get_2d_scalar_shared;
%rename(get_2d_vector) pism::Vars::get_2d_vector_shared;
%rename(get_2d_mask) pism::Vars::get_2d_mask_shared;
%rename(get_2d_cell_type) pism::Vars::get_2d_cell_type_shared;
%rename(get_3d_scalar) pism::Vars::get_3d_scalar_shared;
%rename(keys) pism::Vars::keys_shared;

/* remove Vars::remove */
%ignore pism::Vars::remove;

%extend pism::Vars {
  %pythoncode %{
def get(self, key):
    # try this first (a cell type mask is a mask, which is a scalar)
    try:
        return self.get_2d_cell_type(key)
    except:
        pass

    try:
        return self.get_2d_mask(key)
    except:
        pass

    try:
        return self.get_2d_scalar(key)
    except:
        pass

    try:
        return self.get_2d_vector(key)
    except:
        pass

    try:
        return self.get_3d_scalar(key)
    except:
        pass

    return self._get(key)
  %}
}

%include "util/Vars.hh"
