#
# tlibs2
# @author Tobias Weber <tweber@ill.fr>
# @date 4-jun-2020
# @license see 'LICENSE' file
#

cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

project(tlibs2)


message("Project: ${PROJECT_NAME}, build type: ${CMAKE_BUILD_TYPE}.")

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif()


#set(CMAKE_CXX_STANDARD 17)
#add_definitions(-std=c++17)

set(CMAKE_CXX_STANDARD 20)
add_definitions(-std=c++20)



# -----------------------------------------------------------------------------
# Boost
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED COMPONENTS system filesystem iostreams REQUIRED)
add_definitions(${Boost_CXX_FLAGS})

message("Using Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Python
find_package(PythonLibs)

message("Using Python version ${PYTHONLIBS_VERSION_STRING}.")
# -----------------------------------------------------------------------------



include_directories(
	"${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/libs"
	"${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"${PYTHON_INCLUDE_DIRS}"
)



# -----------------------------------------------------------------------------
# Python module
find_package(SWIG REQUIRED)
set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
include(${SWIG_USE_FILE})

set_source_files_properties(swig_mod/tl2.i PROPERTIES CPLUSPLUS TRUE)
set_source_files_properties(swig_mod/tl2.i PROPERTIES SWIG_FLAGS "-I ${PROJECT_SOURCE_DIR}/libs")

swig_add_library(tl2_py LANGUAGE python SOURCES swig_mod/tl2.i libs/log.cpp)

target_link_libraries(tl2_py "${PYTHON_LIBRARIES}")


# installation
install(TARGETS tl2_py LIBRARY DESTINATION lib64)
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# unit tests
if(BUILD_TESTING)
	add_subdirectory(tests)
endif()
# -----------------------------------------------------------------------------
