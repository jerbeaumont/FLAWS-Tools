set (proj TCLAP)

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location URL https://sourceforge.net/projects/tclap/files/tclap-1.2.4.tar.gz)
endif()

ExternalProject_Add(${proj}
  ${location}
  PREFIX ${CMAKE_BINARY_DIR}/External-Projects/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj}
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  CMAKE_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj}/include)

set(FLAWS_Tools_DEPS "${FLAWS_Tools_DEPS};${proj}")
