cmake_minimum_required(VERSION 2.6)
 
PROJECT(ca_smoothing)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(ca_smoothing ca_smoothing.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(ca_smoothing ${VTK_LIBRARIES})
else()
  target_link_libraries(ca_smoothing vtkHybrid )
endif()
