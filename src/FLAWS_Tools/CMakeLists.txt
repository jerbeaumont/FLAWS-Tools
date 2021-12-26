cmake_minimum_required(VERSION 3.1.0)

cmake_policy(SET CMP0048 NEW)
project(FLAWS_TOOLS VERSION 1.0)

## #############################################################################
## Add additional cmake modules
## #############################################################################

# On Macs, disable MACOSX_RPATH (cmake 3.x)
set(CMAKE_MACOSX_RPATH OFF)

# Use C++ 14 standard
set(CMAKE_CXX_STANDARD 14)

#Add definition for version number
add_definitions(-D${PROJECT_NAME}_VERSION="${${PROJECT_NAME}_VERSION}")

#On windows, add definitions for M_PI
if (WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:strict")
endif()

set(FLAWS_TOOLS_CMAKE_DOCUMENTATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake/documentation)
set(${PROJECT_NAME}_CMAKE_DIRS
  ${CMAKE_SOURCE_DIR}/ # for CTestConfig file
  ${CMAKE_SOURCE_DIR}/cmake
  ${FLAWS_TOOLS_CMAKE_DOCUMENTATION_DIR}
  ${CMAKE_SOURCE_DIR}/cmake/module
  ${CMAKE_CURRENT_BINARY_DIR}/cmake # for auto generated ${PROJECT_NAME}Exports
  )

set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_CMAKE_DIRS}
  ${CMAKE_MODULE_PATH}
  )

include(GenerateExportHeader)
# if GenerateExportHeader is not include explicitly here,
# it may conflict with redefinition of the module in tierce lib (i.e. vtk6)
include(list_source_files)
include(list_header_directories_to_include)
include(set_exe_install_rules)
include(set_lib_install_rules)

################################################################################
# Output directories
################################################################################

set(EXECUTABLE_OUTPUT_PATH
  ${CMAKE_BINARY_DIR}/bin
  CACHE PATH "Single output directory for building all executables."
  )

set(LIBRARY_OUTPUT_PATH
  ${CMAKE_BINARY_DIR}/lib
  CACHE PATH "Single output directory for building all libraries."
  )
  
  
################################################################################
# Setup project
################################################################################

option(BUILD_SHARED_LIBS
  "Build with shared libraries."
  ON
  )

set(${PROJECT_NAME}_LIBRARY_DIRS
  ${LIBRARY_OUTPUT_PATH}
  )

set(${PROJECT_NAME}_USE_FILE
  ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Use.cmake
  )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
configure_file(${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Exports.cmake.in
    ${${PROJECT_NAME}_BINARY_DIR}/cmake/${PROJECT_NAME}Exports.cmake
    @ONLY IMMEDIATE
)

include(${PROJECT_NAME}Exports)
include(FLAWS_ToolsModulesSetup)	

################################################################################
# Find dependencies
################################################################################

include(FLAWS_ToolsDependencies)
include(ITKLibs)

################################################################################
# Setup use file
################################################################################

if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in)
  configure_file( ## Common use file
    ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in
    ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Use.cmake
    @ONLY IMMEDIATE
    )
endif(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in)


################################################################################
# Setup Configuration files
################################################################################

if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in)
  configure_file( ## Build tree configure file
    ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    @ONLY IMMEDIATE
    )
endif(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in)


################################################################################
# add_subdirectories
################################################################################

if (BUILD_MODULE_FLAWS-SEQUENCE)
  add_subdirectory(flaws-sequence)
endif()

if (BUILD_MODULE_SEQUENCE-SIMULATION)
  add_subdirectory(sequence-simulation)
endif()
