cmake_minimum_required(VERSION 2.8)

PROJECT(autoneedle)

FIND_PACKAGE(VTK REQUIRED)

FIND_PACKAGE(GTest REQUIRED)

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
  src/needle.cxx
  src/needlematrix.cxx
  src/sparsediagonalmatrix.cxx
)

INCLUDE_DIRECTORIES(
  gen_src/
  3rdparty/cml-1_0_2/
  include/
  ${GTEST_INCLUDE_DIR}
)

ADD_LIBRARY( needlelib ${source} )
ADD_DEPENDENCIES( needlelib generatedCode )  

ADD_EXECUTABLE( autoneedle src/main.cxx src/Rendering.cxx)
TARGET_LINK_LIBRARIES(autoneedle ${VTK_LIBRARIES})
TARGET_LINK_LIBRARIES(autoneedle needlelib)

ADD_EXECUTABLE( tests tests/testing.cxx)
TARGET_LINK_LIBRARIES( tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} pthread needlelib )
