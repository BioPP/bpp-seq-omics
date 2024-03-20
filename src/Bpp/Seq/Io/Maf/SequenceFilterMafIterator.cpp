// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SequenceFilterMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> SequenceFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_)
  {
    map<string, unsigned int> counts;
    for (size_t i = currentBlock_->getNumberOfSequences(); i > 0; --i)
    {
      string species = currentBlock_->sequence(i - 1).getSpecies();
      if (!VectorTools::contains(species_, species))
      {
        if (logstream_)
        {
          (*logstream_ << "SEQUENCE FILTER: remove sequence '" << species << "' from current block " << currentBlock_->getDescription() << ".").endLine();
        }
        if (!keep_)
        {
          currentBlock_->removeSequence(i - 1);
        }
      }
      else
      {
        counts[species]++;
      }
    }
    bool test = currentBlock_->getNumberOfSequences() == 0;
    // Avoid a memory leak:
    if (test)
    {
      if (logstream_)
      {
        (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " is now empty. Try to get the next one.").endLine();
      }
    }
    else
    {
      test = strict_ && (counts.size() != species_.size());
      if (test)
      {
        if (logstream_)
        {
          (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " does not contain all species and will be ignored. Try to get the next one.").endLine();
        }
      }
      else
      {
        if (rmDuplicates_)
        {
          test = false;
          map<string, unsigned int>::iterator it;
          for (it = counts.begin(); it != counts.end() && !(test = it->second > 1); it++)
          {}
          if (test)
          {
            if (logstream_)
            {
              (*logstream_ << "SEQUENCE FILTER: block " << currentBlock_->getDescription() << " has two sequences for species '" << it->first << "' and will be ignored. Try to get the next one.").endLine();
            }
          }
          else
          {
            return move(currentBlock_);
          }
        }
        else
        {
          return move(currentBlock_);
        }
      }
    }

    // Look for the next block:
    currentBlock_ = iterator_->nextBlock();
  }

  return move(currentBlock_);
}
