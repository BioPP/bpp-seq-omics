// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ChromosomeMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

std::unique_ptr<MafBlock> ChromosomeMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_)
  {
    bool foundRef = false;
    string chr = "";
    for (size_t i = 0; i < currentBlock_->getNumberOfSequences() && !foundRef; ++i)
    {
      string species = currentBlock_->sequence(i).getSpecies();
      if (species == ref_)
      {
        foundRef = true;
        chr = currentBlock_->sequence(i).getChromosome();
      }
    }
    if (!foundRef)
    {
      if (logstream_)
      {
        (*logstream_ << "CHROMOSOME FILTER: block does not contain reference species and was removed.").endLine();
      }
    }
    else if (chr_.find(chr) == chr_.end())
    {
      if (logstream_)
      {
        (*logstream_ << "CHROMOSOME FILTER: reference species without queried chromosome was removed.").endLine();
      }
    }
    else
    {
      return std::move(currentBlock_);
    }

    // Look for the next block:
    currentBlock_ = iterator_->nextBlock();
  }

  return std::move(currentBlock_);
}
