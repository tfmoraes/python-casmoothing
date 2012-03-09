cmake_minimum_required(VERSION 2.6)
 
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")
endif()

PROJECT(ca_smoothing)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(ca_smoothing SHARED ca_smoothing.cxx)
 
add_executable(ca_smoothing_app ca_smoothing_app.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(ca_smoothing_app ca_smoothing ${VTK_LIBRARIES})
else()
  target_link_libraries(ca_smoothing vtkHybrid )
endif()
