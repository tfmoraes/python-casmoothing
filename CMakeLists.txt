cmake_minimum_required(VERSION 2.6)
INCLUDE(CPack)
INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)

PROJECT(ca_smoothing)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
 
find_package(VTK REQUIRED)
find_package(VTK COMPONENTS vtkCommonCore vtkWrappingTools vtkCommonDataModel )
include(${VTK_USE_FILE})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

message(${VTK_LIBRARIES})

#add_library(ca_smoothing SHARED ca_smoothing.cxx)
#add_executable(ca_smoothing_app ca_smoothing_app.cxx)

#SET (Libraries vtkPythonCore)

SET_SOURCE_FILES_PROPERTIES(ca_smoothing.i PROPERTIES CPLUSPLUS ON)
#SET_SOURCE_FILES_PROPERTIES(ca_smoothing.i PROPERTIES SWIG_FLAGS "-includeall")
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
SWIG_ADD_MODULE(ca_smoothing python ca_smoothing.i ca_smoothing.cxx)

message(STATUS "Libraries" ${Libraries})

if(VTK_LIBRARIES)
	#target_link_libraries(ca_smoothing ${VTK_LIBRARIES} ${Libraries})
	#target_link_libraries(ca_smoothing_app ca_smoothing ${VTK_LIBRARIES} ${Libraries})
        SWIG_LINK_LIBRARIES(ca_smoothing ${PYTHON_LIBRARIES} ${Libraries} ${VTK_LIBRARIES})
else()
        #target_link_libraries(ca_smoothing vtkHybrid )
        #target_link_libraries(ca_smoothing_app ca_smoothing vtkHybrid ${Libraries})
	SWIG_LINK_LIBRARIES(ca_smoothing ${PYTHON_LIBRARIES} vtkHybrid vtkCommon vtkHybridPythonD vtkFiltering)
       
endif()

IF (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

	# If g++ is greater than 4.6, set to use c++0x standard.
	execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GXX_VERSION)
	if(GXX_VERSION VERSION_GREATER "4.6")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	endif(GXX_VERSION VERSION_GREATER "4.6")

ENDIF(CMAKE_COMPILER_IS_GNUCXX)

CHECK_INCLUDE_FILE_CXX (unordered_map HAVE_CXX0X)
CHECK_INCLUDE_FILE_CXX (tr1/unordered_map HAVE_TR1)

IF (NOT HAVE_CXX0X)
	SET (HAVE_CXX0X 0)
ENDIF (NOT HAVE_CXX0X)

IF (NOT HAVE_TR1)
	SET (HAVE_TR1 0)
ENDIF (NOT HAVE_TR1)


# Finding python site-packages dir
execute_process( COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(prefix='')" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE ) 

INSTALL(TARGETS _ca_smoothing DESTINATION ${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_PACKAGES})
INSTALL(FILES ${CMAKE_BINARY_DIR}/lib/ca_smoothing.py DESTINATION ${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_PACKAGES})
