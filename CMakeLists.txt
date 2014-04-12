cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(VTK REQUIRED)

INCLUDE(${VTK_USE_FILE})


ADD_CUSTOM_COMMAND(
  OUTPUT
  ${PROJECT_SOURCE_DIR}/gen_src/generatedCode.h
  COMMAND
  python ${PROJECT_SOURCE_DIR}/tools/GenerateCode.py
  WORKING_DIRECTORY
  ${PROJECT_SOURCE_DIR}/
  DEPENDS
  ${PROJECT_SOURCE_DIR}/tools/GenerateCode.py
  ${PROJECT_SOURCE_DIR}/tools/CodeToC.py
)

ADD_CUSTOM_TARGET(
  generatedCode
  ALL
  DEPENDS 
  ${PROJECT_SOURCE_DIR}/gen_src/generatedCode.h
)

SET( source
  src/main.cxx
  src/needle.cxx
  src/sparsediagonalmatrix.cxx
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
