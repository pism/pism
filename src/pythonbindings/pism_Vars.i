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

%extend pism::Vars
{
  %pythoncode "Vars.py"
}

%include "util/Vars.hh"
