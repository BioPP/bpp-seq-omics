# CMake script for Bio++ Sequence Omics Library
# Authors:
#   Sylvain Gaillard
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 30/07/2012

# File list
set (CPP_FILES
  Bpp/Seq/Feature/Bed/BedGraphFeatureReader.cpp
  Bpp/Seq/Feature/Gff/GffFeatureReader.cpp
  Bpp/Seq/Feature/Gtf/GtfFeatureReader.cpp
  Bpp/Seq/Feature/SequenceFeature.cpp
  Bpp/Seq/Feature/SequenceFeatureTools.cpp
  Bpp/Seq/Io/Fastq.cpp
  Bpp/Seq/Io/Maf/AlignmentFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/BlockMergerMafIterator.cpp
  Bpp/Seq/Io/Maf/ChromosomeMafIterator.cpp
  Bpp/Seq/Io/Maf/ConcatenateMafIterator.cpp
  Bpp/Seq/Io/Maf/CoordinateTranslatorMafIterator.cpp
  Bpp/Seq/Io/Maf/CoordinatesOutputMafIterator.cpp
  Bpp/Seq/Io/Maf/DuplicateFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/EntropyFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/FeatureExtractorMafIterator.cpp
  Bpp/Seq/Io/Maf/FeatureFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/FullGapFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/IterationListener.cpp
  Bpp/Seq/Io/Maf/MafIterator.cpp
  Bpp/Seq/Io/Maf/MafParser.cpp
  Bpp/Seq/Io/Maf/MafSequence.cpp
  Bpp/Seq/Io/Maf/MafStatistics.cpp
  Bpp/Seq/Io/Maf/MaskFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/MsmcOutputMafIterator.cpp
  Bpp/Seq/Io/Maf/TableOutputMafIterator.cpp
  Bpp/Seq/Io/Maf/OrphanSequenceFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/OutputAlignmentMafIterator.cpp
  Bpp/Seq/Io/Maf/OutputMafIterator.cpp
  Bpp/Seq/Io/Maf/PlinkOutputMafIterator.cpp
  Bpp/Seq/Io/Maf/QualityFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/RemoveEmptySequencesMafIterator.cpp
  Bpp/Seq/Io/Maf/SequenceFilterMafIterator.cpp
  Bpp/Seq/Io/Maf/SequenceLDhotOutputMafIterator.cpp
  Bpp/Seq/Io/Maf/SequenceStatisticsMafIterator.cpp
  Bpp/Seq/Io/Maf/SequenceStreamToMafIterator.cpp
  Bpp/Seq/Io/Maf/VcfOutputMafIterator.cpp
  Bpp/Seq/Io/Maf/WindowSplitMafIterator.cpp
  )

# Build the static lib
add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC})

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED})

# Install libs and headers
install (
  TARGETS ${PROJECT_NAME}-static ${PROJECT_NAME}-shared
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
install (
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Bpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
  )
# Install cmake file representing targets
install (EXPORT ${PROJECT_NAME}-targets DESTINATION ${cmake-package-location})
