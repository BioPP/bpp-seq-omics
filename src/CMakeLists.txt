# CMake script for Bio++ Sequence Library
# Author: Sylvain Gaillard and Julien Dutheil
# Created: 30/07/2012

# File list
set(CPP_FILES
  Bpp/Seq/Io/Fastq.cpp
  Bpp/Seq/Io/Maf/MafAlignmentParser.cpp
  Bpp/Seq/Io/Maf/MafIterator.cpp
  Bpp/Seq/Io/Maf/MafSequence.cpp
  Bpp/Seq/Io/Maf/MafStatistics.cpp
  Bpp/Seq/Io/Maf/IterationListener.cpp
  Bpp/Seq/Feature/SequenceFeature.cpp
  Bpp/Seq/Feature/Gff/GffFeatureReader.cpp
  Bpp/Seq/Feature/Gtf/GtfFeatureReader.cpp
  )

set(H_FILES
  Bpp/Seq/Io/Fastq.h
  Bpp/Seq/Io/Maf/MafSequence.h
  Bpp/Seq/Io/Maf/MafAlignmentParser.h
  Bpp/Seq/Io/Maf/MafIterator.h
  Bpp/Seq/Io/Maf/MafStatistics.h
  Bpp/Seq/Io/Maf/IterationListener.h
  Bpp/Seq/Feature/SequenceFeature.h
  Bpp/Seq/Feature/Gff/GffFeatureReader.h
  Bpp/Seq/Feature/Gtf/GtfFeatureReader.h
  )

# Build the static lib
add_library(bppseqomics-static STATIC ${CPP_FILES})
set_target_properties(bppseqomics-static
  PROPERTIES OUTPUT_NAME bpp-seq-omics
  CLEAN_DIRECT_OUTPUT 1
  )
target_link_libraries(bppseqomics-static ${LIBS})

# Build the shared lib
add_library(bppseqomics-shared SHARED ${CPP_FILES})
set_target_properties(bppseqomics-shared
  PROPERTIES OUTPUT_NAME bpp-seq-omics
  CLEAN_DIRECT_OUTPUT 1
  VERSION ${BPPSEQ_VERSION}
  SOVERSION ${BPPSEQ_VERSION_MAJOR}
  )
target_link_libraries(bppseqomics-shared ${LIBS})

# Install libs
install(TARGETS bppseqomics-static bppseqomics-shared DESTINATION lib${LIB_SUFFIX})

# Install headers
install(DIRECTORY Bpp/ DESTINATION include/Bpp FILES_MATCHING PATTERN "*.h" PATTERN "*.all")
