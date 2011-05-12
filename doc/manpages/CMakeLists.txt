find_program(RST2MAN NAMES rst2man rst2man.py rst2man-2.6.py)
mark_as_advanced(RST2MAN)

if (RST2MAN)

  macro(add_manpages names)
    foreach(name ${names})
      add_custom_command(OUTPUT ${name}.1
        COMMAND ${RST2MAN} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.rst ${name}.1
        DEPENDS ${name}.rst
        WORKING_DIRECTORY)

      install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${name}.1
        DESTINATION "share/man/man1")
    endforeach(name)
  endmacro(add_manpages)

  add_manpages("pismr;pisms;pismv;pgrn;pross")

  add_custom_target(manpages
    DEPENDS pismr.1 pisms.1 pismv.1 pgrn.1 pross.1)

endif()