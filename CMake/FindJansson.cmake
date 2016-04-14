# Copyright (C) 2014 Johannes Schauer <j.schauer@email.de>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

find_path(JANSSON_INCLUDE_DIRS NAMES "jansson.h")
find_library(JANSSON_LIBRARIES NAMES "jansson")
mark_as_advanced(JANSSON_INCLUDE_DIRS JANSSON_LIBRARIES)

if (JANSSON_INCLUDE_DIRS)
	file(STRINGS "${JANSSON_INCLUDE_DIRS}/jansson.h" jansson_version_str REGEX "^#define[\t ]+JANSSON_VERSION[\t ]+\".*\"")
	string(REGEX REPLACE "^#define[\t ]+JANSSON_VERSION[\t ]+\"([^\"]*)\".*" "\\1" JANSSON_VERSION_STRING "${jansson_version_str}")
	unset(jansson_version_str)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Jansson
	REQUIRED_VARS JANSSON_LIBRARIES JANSSON_INCLUDE_DIRS
	VERSION_VAR JANSSON_VERSION_STRING)
