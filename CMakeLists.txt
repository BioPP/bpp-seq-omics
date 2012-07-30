# CMake script for bpp-seq unit tests
# Author: Julien Dutheil
# Created: 30/07/2012

MACRO(TEST_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})

  SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    MESSAGE(WARNING "Library ${lib_name} is already installed in the system tree. Test will be built against it. This may lead to unexpected results. You may want to do 'make install' before 'make test', or remove the installed version.")
  ELSE()
    SET(${lib_name}_LIBRARY "-L../src -lbpp-seq")
    SET(${lib_name}_INCLUDE_DIR "../src/")
  ENDIF()
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(TEST_FIND_LIBRARY)

#Find the bpp-seq library library:
TEST_FIND_LIBRARY(LIBS bpp-seq-omics Bpp/Seq/Feature/SequenceFeature.h)

ADD_EXECUTABLE(test_gff test_gff.cpp)
TARGET_LINK_LIBRARIES(test_gff ${LIBS})
ADD_TEST(test_gff "test_gff")

ADD_EXECUTABLE(test_fastq test_fastq.cpp)
TARGET_LINK_LIBRARIES(test_fastq ${LIBS})
ADD_TEST(test_fastq "test_fastq")

IF(UNIX)
  SET_PROPERTY(TEST test_gff test_fastq PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:../src")
ENDIF()

IF(APPLE)
  SET_PROPERTY(TEST test_gff test_fastq PROPERTY ENVIRONMENT "DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}:../src")
ENDIF()

IF(WIN32)
  SET(ENV{PATH} "$ENV{PATH};..\\src")
ENDIF()

