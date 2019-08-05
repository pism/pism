
find_program(CUT_EXECUTABLE cut)
find_program(XARGS_EXECUTABLE xargs)
mark_as_advanced(CUT_EXECUTABLE XARGS_EXECUTABLE)

add_subdirectory(debian)
