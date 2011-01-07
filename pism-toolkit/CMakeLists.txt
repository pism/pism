find_program (FAKEROOT fakeroot)
find_program (DEBHELPER dh)

mark_as_advanced (FORCE FAKEROOT)
mark_as_advanced (FORCE DEBHELPER)

set (toolkit_version "0.1")
set (toolkit_dir pism-toolkit-${toolkit_version})

if (FAKEROOT MATCHES "FAKEROOT-NOTFOUND" OR DEBHELPER MATCHES "DEBHELPER-NOTFOUND")
  message (STATUS "Either fakeroot or debhelper was not found. Debian package targets are disabled.")
else ()
  foreach (file changelog compat control copyright rules)
    configure_file (debian/${file} ${toolkit_dir}/debian/${file} COPYONLY)
  endforeach ()

  add_custom_target (pism-toolkit
    COMMAND fakeroot debian/rules binary
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${toolkit_dir})
endif ()
