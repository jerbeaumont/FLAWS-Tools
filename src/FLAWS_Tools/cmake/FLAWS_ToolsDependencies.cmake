# ITK
find_package(ITK REQUIRED)

set(_ITKVersionString "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}" )
math(EXPR _ITKVersionNum "${ITK_VERSION_MAJOR}*100*100 + ${ITK_VERSION_MINOR}*100 + ${ITK_VERSION_PATCH}")
  
if(_ITKVersionNum LESS 41001)
  message(SEND_ERROR "The ITK version you want to use (${_ITKVersionString}) is not supported by this project. Please use a more recent version of ITK. The minimum required version is 4.10.1")
else()
  include(${ITK_USE_FILE})
endif()

# TCLAP
option(BUILD_TOOLS "Build command line executables" ON)
mark_as_advanced(BUILD_TOOLS)
if(BUILD_TOOLS)
  find_package(TCLAP REQUIRED)
endif()
