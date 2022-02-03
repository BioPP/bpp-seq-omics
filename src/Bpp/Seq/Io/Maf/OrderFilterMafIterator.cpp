//
// File: OrderFilterMafIterator.cpp
// Authors: Julien Dutheil
// Created: Sat Jan 20 2018
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2018)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

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
  if (!block.hasMafSequenceForSpecies(refSpecies_))
    return true; // We consider a block with no reference sequence as ordered
  const MafSequence& refSeq = block.getMafSequenceForSpecies(refSpecies_);
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
