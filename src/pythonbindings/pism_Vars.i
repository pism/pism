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
