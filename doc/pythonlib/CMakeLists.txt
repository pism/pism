set(SPHYNX_BUILD_EXECUTABLE sphinx-build)

set(RST_FILES
index.rst 
api/index.rst
api/PISM.rst
api/ssa.rst
api/sia.rst
api/invert.rst
api/invert_ssa.rst
api/invert_ssa_tao.rst
api/invert_ssa_siple.rst
api/sipletools.rst

)

set (sphinx_helpers
_static/pymunk.js
pismdoc/theme.conf
pismdoc/static/PISM.css )


# copy all the files needed to generate browser sources
foreach (filename ${sphinx_helpers})
  configure_file(${filename} ${filename} COPYONLY)
endforeach()

foreach (filename ${RST_FILES})
  configure_file(${filename} ${filename} COPYONLY)
endforeach()

configure_file(conf.py conf.py @ONLY)


add_custom_target (sphinx_docs
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/doc/pythonlib/_static ${CMAKE_BINARY_DIR}/doc/pythonlib/_static
  COMMAND ${SPHYNX_BUILD_EXECUTABLE} -b html . html
  DEPENDS ${RST_FILES} ${sphinx_helpers} install_python_bindings_in_the_build_directory
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
