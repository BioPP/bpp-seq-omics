// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "DuplicateFilterMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> DuplicateFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_)
  {
    bool foundRef = false;
    string chr = "";
    char strand = '+';
    size_t start = 0;
    size_t stop  = 0;
    for (size_t i = 0; i < currentBlock_->getNumberOfSequences() && !foundRef; ++i)
    {
      string species = currentBlock_->sequence(i).getSpecies();
      if (species == ref_)
      {
        foundRef = true;
        chr    = currentBlock_->sequence(i).getChromosome();
        strand = currentBlock_->sequence(i).getStrand();
        start  = currentBlock_->sequence(i).start();
        stop   = currentBlock_->sequence(i).stop();
      }
    }
    if (!foundRef)
    {
      if (logstream_)
      {
        (*logstream_ << "DUPLICATE FILTER: block does not contain reference species and was removed.").endLine();
      }
    }
    else
    {
      size_t occurrence = blocks_[chr][strand][start][stop]++;
      if (occurrence > 0)
      {
        if (logstream_)
        {
          (*logstream_ << "DUPLICATE FILTER: sequence in reference species was found in a previous block. New block was removed.").endLine();
        }
      }
      else
      {
        return std::move(currentBlock_);
      }
    }

    // Look for the next block:
    currentBlock_ = iterator_->nextBlock();
  }

  return std::move(currentBlock_);
}
