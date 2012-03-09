cmake_minimum_required(VERSION 2.6)
 
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")
endif()

PROJECT(ca_smoothing)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})


add_library(ca_smoothing SHARED ca_smoothing.cxx)
add_executable(ca_smoothing_app ca_smoothing_app.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(ca_smoothing ${VTK_LIBRARIES})
  target_link_libraries(ca_smoothing_app ca_smoothing ${VTK_LIBRARIES})
else()
  target_link_libraries(ca_smoothing vtkHybrid )
endif()

SET_SOURCE_FILES_PROPERTIES(ca_smoothing.i PROPERTIES CPLUSPLUS ON)
#SET_SOURCE_FILES_PROPERTIES(ca_smoothing.i PROPERTIES SWIG_FLAGS "-includeall")

SWIG_ADD_MODULE(ca_smoothing python ca_smoothing.i ca_smoothing.cxx)
SWIG_LINK_LIBRARIES(ca_smoothing ${PYTHON_LIBRARIES} ${VTK_LIBRARIES})
