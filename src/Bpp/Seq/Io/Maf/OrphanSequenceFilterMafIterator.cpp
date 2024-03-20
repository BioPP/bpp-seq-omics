// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OrphanSequenceFilterMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

std::unique_ptr<MafBlock> OrphanSequenceFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_)
  {
    map<string, unsigned int> counts;
    for (size_t i = 0; i < currentBlock_->getNumberOfSequences(); ++i)
    {
      string species = currentBlock_->sequence(i).getSpecies();
      counts[species]++;
    }
    bool test = counts.size() <= species_.size();
    if (test)
    {
      // We have to check that the species are the right one:
      bool loseCrit = false;
      bool strictCrit = true;
      bool duplicate = false;
      for (size_t i = 0; i < species_.size() && !duplicate; ++i)
      {
        map<string, unsigned int>::iterator it = counts.find(species_[i]);
        if (it != counts.end())
        {
          loseCrit = true;
          if (rmDuplicates_ && it->second > 1)
          {
            // Duplicated sequences, block is discarded if asked to...
            duplicate = true;
          }
        }
        else
        {
          strictCrit = false;
        }
      }
      if (!duplicate)  // No duplicate found or duplicates are kept
      {
        if (strictCrit || (!strict_ && loseCrit))
          return std::move(currentBlock_);
      }
    }
    // Otherwise there is at least one extra species, we get the next block...
    currentBlock_ = iterator_->nextBlock();
  }

  return std::move(currentBlock_);
}
