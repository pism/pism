# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

mark_as_advanced(CGAL_DIR CGAL_DONT_OVERRIDE_CMAKE_FLAGS)

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

add_to_cached_list( CGAL_EXECUTABLE_TARGETS naturalNeighbor )

include_directories(${FEvoR_INCLUDES})

add_library(pismfevor
  PISMFEvoR.cc
  FEvoR_IO.cc)

target_link_libraries(pismfevor pismutil pismstressbalance ${FEvoR_LIB} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismfevor DESTINATION ${Pism_LIB_DIR})
endif()
