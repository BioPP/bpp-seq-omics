// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OrderFilterMafIterator.h"

// From bpp-seq:
// #include <Bpp/Seq/SequenceWithAnnotationTools.h>
// #include <Bpp/Seq/Container/VectorSiteContainer.h>
// #include <Bpp/Seq/SiteTools.h>
// #include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

bool OrderFilterMafIterator::parseBlock_(const MafBlock& block)
{
  // Get the reference species for coordinates:
  if (!block.hasSequenceForSpecies(refSpecies_))
    return true; // We consider a block with no reference sequence as ordered
  const auto& refSeq = block.sequenceForSpecies(refSpecies_);
  string chr = refSeq.getChromosome();
  if (chr != currentChr_)
  {
    currentChr_ = chr;
    previousBlockStart_ = 0;
    previousBlockStop_ = 0;
  }
  else
  {
    // Check that block are ordered according to reference sequence:
    if (refSeq.start() < previousBlockStop_)
    {
      if (refSeq.stop() < previousBlockStart_)
      {
        // non-overlapping
        if (unsortedBlockThrowsException_)
          throw Exception("OrderFilterMafIterator: blocks are not ordered according to reference sequence: " + refSeq.getDescription() + "<!>" + TextTools::toString(previousBlockStart_) + ".");
        if (unsortedBlockDiscarded_)
        {
          if (logstream_)
          {
            (*logstream_ << "ORDER FILTER: block " << refSeq.getDescription() << " is not sorted according to previous block and was discarded. Previous block started at " << previousBlockStart_ << " and ended at " << previousBlockStop_ << ".").endLine();
          }
          return false;
        }
      }
      else
      {
        // overlapping
        if (overlappingBlockThrowsException_)
          throw Exception("OrderFilterMafIterator: blocks are overlapping according to reference sequence: " + refSeq.getDescription() + "<!>" + TextTools::toString(previousBlockStop_) + ".");
        if (overlappingBlockDiscarded_)
        {
          if (logstream_)
          {
            (*logstream_ << "ORDER FILTER: block " << refSeq.getDescription() << " is overlapping with previous block and was discarded. Previous block started at " << previousBlockStart_ << " and ended at " << previousBlockStop_ << ".").endLine();
          }
          return false;
        }
      }
    }
    previousBlockStart_ = refSeq.start();
    previousBlockStop_ = refSeq.stop();
  }
  return true;
}
