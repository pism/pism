
cmake_minimum_required(VERSION 2.8)

PROJECT(WriteVTU)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


add_executable(WriteVTU MACOSX_BUNDLE WriteVTU)

if(VTK_LIBRARIES)
  target_link_libraries(WriteVTU ${VTK_LIBRARIES})
else()
  target_link_libraries(WriteVTU vtkHybrid vtkWidgets)
endif()

target_link_libraries (WriteVTU netcdf_c++4)
