# CMake script for Bio++ Sequence Library
# Author: Sylvain Gaillard and Julien Dutheil
# Created: 30/07/2012

# File list
SET(CPP_FILES
  Bpp/Seq/Io/Fastq.cpp
  Bpp/Seq/Io/Maf/MafAlignmentParser.cpp
  Bpp/Seq/Io/Maf/MafIterator.cpp
  Bpp/Seq/Io/Maf/MafSequence.cpp
  Bpp/Seq/Io/Maf/MafStatistics.cpp
  Bpp/Seq/Io/Maf/IterationListener.cpp
  Bpp/Seq/Feature/SequenceFeature.cpp
  Bpp/Seq/Feature/SequenceFeatureTools.cpp
  Bpp/Seq/Feature/Gff/GffFeatureReader.cpp
  Bpp/Seq/Feature/Gtf/GtfFeatureReader.cpp
  )

SET(H_FILES
  Bpp/Seq/Io/Fastq.h
  Bpp/Seq/Io/Maf/MafSequence.h
  Bpp/Seq/Io/Maf/MafAlignmentParser.h
  Bpp/Seq/Io/Maf/MafIterator.h
  Bpp/Seq/Io/Maf/MafStatistics.h
  Bpp/Seq/Io/Maf/IterationListener.h
  Bpp/Seq/Feature/SequenceFeature.h
  Bpp/Seq/Feature/SequenceFeatureTools.h
  Bpp/Seq/Feature/Gff/GffFeatureReader.h
  Bpp/Seq/Feature/Gtf/GtfFeatureReader.h
  )

# Build the static lib
ADD_LIBRARY(bppseqomics-static STATIC ${CPP_FILES})
SET_TARGET_PROPERTIES(bppseqomics-static
  PROPERTIES OUTPUT_NAME bpp-seq-omics
  CLEAN_DIRECT_OUTPUT 1
  )
TARGET_LINK_LIBRARIES(bppseqomics-static ${LIBS})

# Build the shared lib
ADD_LIBRARY(bppseqomics-shared SHARED ${CPP_FILES})
SET_TARGET_PROPERTIES(bppseqomics-shared
  PROPERTIES OUTPUT_NAME bpp-seq-omics
  CLEAN_DIRECT_OUTPUT 1
  VERSION ${BPPSEQOMICS_VERSION}
  SOVERSION ${BPPSEQOMICS_VERSION_MAJOR}
  )
TARGET_LINK_LIBRARIES(bppseqomics-shared ${LIBS})

# Install libs
INSTALL(TARGETS bppseqomics-static bppseqomics-shared DESTINATION lib${LIB_SUFFIX})

# Install headers
INSTALL(DIRECTORY Bpp/ DESTINATION include/Bpp FILES_MATCHING PATTERN "*.h" PATTERN "*.all")

INSTALL(CODE "EXECUTE_PROCESS(COMMAND python ${CMAKE_SOURCE_DIR}/genIncludes.py ${CMAKE_PREFIX_PATH}/include/Bpp/Seq)")

