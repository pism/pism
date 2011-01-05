find_program (FAKEROOT fakeroot)
find_program (DPKG dpkg)

if (FAKEROOT MATCHES "FAKEROOT-NOTFOUND" OR DPKG MATCHES "DPKG-NOTFOUND")
  message (STATUS "Either fakeroot or dpkg was not found. Debian package targets are disabled.")
else ()
  # If we're in a Subversion working copy, export to a temporary directory:
  if (EXISTS ${Pism_SOURCE_DIR}/.svn)
    add_custom_target (pism-toolkit
      COMMAND svn export ${Pism_SOURCE_DIR}/debian/pism-toolkit ${PROJECT_BINARY_DIR}/pism-toolkit
      COMMAND fakeroot dpkg --build pism-toolkit
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )
  else ()
    add_custom_target (pism-toolkit
      COMMAND fakeroot dpkg --build ${Pism_SOURCE_DIR}/debian/pism-toolkit
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )
  endif ()
endif ()
