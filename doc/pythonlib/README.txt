A temporary readme for now.

Documentation building should be added to the CMAKE build.  For now:

1) Install sphinx (http://sphinx-doc.org)
2) Build PISM with python bindings.
3) Add PISM to your PYTHONPATH
4) sphinx-build -b html [this directory] build/html

   where [this directory]  is $Pism_SOURCE_DIR/doc/pythonlib.

   Ignore the warnings:

   WARNING! There are options you set that were not used!
   WARNING! could be spelling mistake, etc!
   Option left: name:-b value: html
   
   which are odd side effects of sphinx invoking python 
  invoking PISM including PETSc and PETSc complaining 
  about unused command line options.

Documentation root will be built in build/html/index.html

