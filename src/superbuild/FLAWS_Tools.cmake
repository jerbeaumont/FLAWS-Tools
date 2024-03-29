set (proj FLAWS_Tools)

set (cmake_args
  ${common_cache_args}
  -DEXECUTABLE_OUTPUT_PATH=${CMAKE_BINARY_DIR}/bin
  -DLIBRARY_OUTPUT_PATH=${CMAKE_BINARY_DIR}/lib
  -DBUILD_TOOLS:BOOL=${BUILD_FLAWS_TOOLS_TOOLS}
  -DBUILD_ALL_MODULES:BOOL=ON
  -DBUILD_MODULE_FLAWS-SEQUENCE:BOOL=ON
  -DBUILD_MODULE_SEQUENCE-SIMULATION:BOOL=ON
  -DITK_DIR:PATH=${ITK_BUILD_DIR}
  -DTCLAP_INCLUDE_DIR:PATH=${TCLAP_SRC_DIR}
  )

ExternalProject_Add(${proj}
  DEPENDS ${${proj}_DEPS}
  PREFIX ${CMAKE_BINARY_DIR}/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/${proj}
  CMAKE_GENERATOR ${cmake_gen}
  CMAKE_ARGS ${cmake_args}
  BUILD_ALWAYS 1
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(${proj} binary_dir)
set(${proj}_BUILD_DIR ${binary_dir})
