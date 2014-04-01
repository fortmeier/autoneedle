cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(VTK REQUIRED)

INCLUDE(${VTK_USE_FILE})

ADD_CUSTOM_TARGET(
  generatedCode
  ALL
  COMMAND python tools/GenerateCode.py
  DEPENDS 
  tools/GenerateCode.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

SET( source
  src/main.cxx
  src/Rendering.cxx
)

INCLUDE_DIRECTORIES(
  gen_src/
  3rdparty/cml-1_0_2/
  include/
)

ADD_EXECUTABLE( autoneedle ${source} )
ADD_DEPENDENCIES( autoneedle generatedCode )  
TARGET_LINK_LIBRARIES(autoneedle ${VTK_LIBRARIES})
