#
# tlibs2 -- tests
# @author Tobias Weber <tweber@ill.fr>
# @date 8-jun-2020
# @license see 'LICENSE' file
#

cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)
enable_testing()

project(tlibs2_tests)


message("Build type: ${CMAKE_BUILD_TYPE}")

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif()


set(CMAKE_CXX_STANDARD 20)
add_definitions(-std=c++20)


# -----------------------------------------------------------------------------
# Boost
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED)
add_definitions(${Boost_CXX_FLAGS})

message("Using Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")
# -----------------------------------------------------------------------------


include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/..")


# -----------------------------------------------------------------------------
# build individual test cases
add_executable(cryst cryst.cpp)
add_executable(mat1 mat1.cpp)
add_executable(stat1 stat1.cpp)
add_executable(quat1 quat1.cpp)
add_executable(leastsq leastsq.cpp)
add_executable(expr expr.cpp ../libs/log.cpp)
#add_executable(eig eig.cpp)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# run individual test cases
add_test(cryst cryst)
add_test(mat1 mat1)
add_test(stat1 stat1)
add_test(quat1 quat1)
add_test(lestsq leastsq)
add_test(expr expr)
#add_test(eig eig)
# -----------------------------------------------------------------------------
