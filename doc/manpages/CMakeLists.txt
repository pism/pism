find_program(RST2MAN_EXECUTABLE NAMES rst2man rst2man.py)
mark_as_advanced(RST2MAN_EXECUTABLE)

if (RST2MAN_EXECUTABLE)
  find_package(UnixCommands)
  macro(add_manpages names)
    foreach(name ${names})
      # Note the sed command: this fixes formatting of command-line options: rst2man
      # cannot handle PISM's long options with one dash.
      add_custom_command(OUTPUT ${name}.1.gz
        COMMAND ${RST2MAN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.rst ${name}.1
        COMMAND sed -i -E -e "s/(-[a-zA-Z]+)\\\\fB /\\1/" ${name}.1
        COMMAND ${GZIP} -n -9 ${name}.1
        DEPENDS ${name}.rst
        VERBATIM)

      install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${name}.1.gz
        DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
    endforeach(name)
  endmacro(add_manpages)

  add_manpages("pismr")

  add_custom_target(manpages ALL DEPENDS pismr.1.gz)
endif()
