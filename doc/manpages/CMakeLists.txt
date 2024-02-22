find_program(RST2MAN_EXECUTABLE NAMES rst2man rst2man.py)
mark_as_advanced(RST2MAN_EXECUTABLE)

if (RST2MAN_EXECUTABLE)
  macro(add_manpages names)
    foreach(name ${names})
      # Note the sed command: this fixes formatting of command-line options: rst2man
      # cannot handle PISM's long options with one dash.
      add_custom_command(OUTPUT ${name}.1
        COMMAND ${RST2MAN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.rst ${name}.1
        COMMAND sed -i -E -e "s/(-[a-zA-Z]+)\\\\fB /\\1/" ${name}.1
        DEPENDS ${name}.rst
        WORKING_DIRECTORY
        VERBATIM)

      install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${name}.1
        DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
    endforeach(name)
  endmacro(add_manpages)

  add_manpages("pismr;pismv")

  add_custom_target(manpages DEPENDS pismr.1 pismv.1)
endif()
