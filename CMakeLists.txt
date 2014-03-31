cmake_minimum_required(VERSION 2.8)

#ADD_CUSTOM_COMMAND(
#  OUTPUT ${CMAKE_SOURCE_DIR}/gen_src/generatedCode.h
#  COMMAND python tools/GenerateCode.py
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/
#)


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
)

INCLUDE_DIRECTORIES(
  gen_src
)

ADD_EXECUTABLE( autoneedle ${source} )
  
