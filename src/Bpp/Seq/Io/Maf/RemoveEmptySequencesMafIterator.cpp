// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "RemoveEmptySequencesMafIterator.h"

using namespace bpp;
using namespace std;

unique_ptr<MafBlock> RemoveEmptySequencesMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_)
  {
    for (size_t i = currentBlock_->getNumberOfSequences(); i > 0; --i)
    {
      const MafSequence& seq = currentBlock_->sequence(i - 1);
      bool isEmpty = true;
      if (unresolvedAsGaps_)
      {
        for (size_t j = 0; isEmpty && j < currentBlock_->getNumberOfSites(); ++j)
        {
          if (!AlphabetTools::DNA_ALPHABET->isUnresolved(seq[j]) && !AlphabetTools::DNA_ALPHABET->isGap(seq[j]))
            isEmpty = false;
        }
      }
      else
      {
        for (size_t j = 0; isEmpty && j < currentBlock_->getNumberOfSites(); ++j)
        {
          if (!AlphabetTools::DNA_ALPHABET->isGap(seq[j]))
            isEmpty = false;
        }
      }
      if (isEmpty)
      {
        currentBlock_->removeSequence(i - 1);
      }
    }
  }
  return std::move(currentBlock_);
}
