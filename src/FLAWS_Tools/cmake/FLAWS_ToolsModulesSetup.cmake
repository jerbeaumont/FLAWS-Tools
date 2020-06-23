OPTION(BUILD_ALL_MODULES "Build all FLAWS_TOOLS modules" ON)

OPTION(BUILD_MODULE_FLAWS-SEQUENCE "Build FLAWS_TOOLS flaws-sequence module" ON)

# Cascading cases to make sure the right modules are indeed turned on

if (BUILD_ALL_MODULES)
  set(BUILD_MODULE_FLAWS-SEQUENCE ON CACHE BOOL "Build FLAWS_TOOLS flaws-sequence module" FORCE)
endif()

set (${PROJECT_NAME}_HEADER_PATHS "")

set(${PROJECT_NAME}_INCLUDE_DIRS 
  "${TCLAP_INCLUDE_DIR}"
)

if (BUILD_MODULE_FLAWS-SEQUENCE)
  list(APPEND ${PROJECT_NAME}_HEADER_PATHS ${${PROJECT_NAME}_SOURCE_DIR}/flaws-sequence)
endif()

list_header_directories_to_include(${PROJECT_NAME} ${${PROJECT_NAME}_HEADER_PATHS})
include_directories(${${PROJECT_NAME}_INCLUDE_DIRS})
