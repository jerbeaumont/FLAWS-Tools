cmake_minimum_required(VERSION 3.1.0)

cmake_policy(SET CMP0048 NEW)
project(FLAWS_TOOLS_SUPERBUILD VERSION 1.0)

set(${PROJECT_NAME}_CMAKE_DIRS
   ${CMAKE_SOURCE_DIR}/FLAWS_Tools/cmake
   ${CMAKE_SOURCE_DIR}/superbuild
   )

set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_CMAKE_DIRS}
  ${CMAKE_MODULE_PATH}
  )

set (CMAKE_CXX_STANDARD 14)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_FLAWS_TOOLS_TOOLS "Build FLAWS_TOOLS tools" ON)

# Use github SSH connection or https ?

option(USE_GITHUB_SSH 
  "Use by default Git SSH addresses, requires public key set on github" OFF
  )
mark_as_advanced(USE_GITHUB_SSH)

include(CheckEnvironment)
# Deal with external dependencies

include(ExternalProject)
include(ConfigureExternalProjectVariables)

#Variables for dependencies and update custom, updated in external projects
set(FLAWS_Tools_DEPS "")
set(Update_Repositories "")

# TCLAP
option(USE_SYSTEM_TCLAP "Use system installed TCLAP" OFF)

if (USE_SYSTEM_TCLAP)
  find_package(TCLAP REQUIRED)
  set(TCLAP_SRC_DIR ${TCLAP_INCLUDE_DIR})
else()
  include(TCLAP)
endif()

# ITK
option(USE_SYSTEM_ITK "Use system installed ITK" OFF)

if (USE_SYSTEM_ITK)
  find_package(ITK REQUIRED)
  set(ITK_BUILD_DIR ${ITK_DIR})
else()
  include(ITK)
endif()

include(FLAWS_Tools)

configure_file(${CMAKE_SOURCE_DIR}/superbuild/BinariesPackaging.cmake.in
  ${CMAKE_BINARY_DIR}/BinariesPackaging.cmake
  @ONLY IMMEDIATE)

add_custom_target(pack-binaries
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/BinariesPackaging.cmake
  COMMENT "-- Generating binary packages"
  )
