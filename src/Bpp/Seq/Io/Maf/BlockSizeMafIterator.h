// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _BLOCKSIZEMAFITERATOR_H_
#define _BLOCKSIZEMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to keep only the ones with a minimum number of species.
 */
class BlockSizeMafIterator :
  public AbstractFilterMafIterator
{
private:
  unsigned int minSize_;

public:
  BlockSizeMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      unsigned int minSize) :
    AbstractFilterMafIterator(iterator),
    minSize_(minSize)
  {}

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_() override
  {
    bool test;
    do
    {
      currentBlock_ = iterator_->nextBlock();
      if (!currentBlock_) break;
      test = (currentBlock_->getNumberOfSequences() < minSize_);
      if (test)
      {
        if (logstream_)
        {
          (*logstream_ << "BLOCK SIZE FILTER: block " << currentBlock_->getDescription() << " with size " << currentBlock_->getNumberOfSites() << " was discarded.").endLine();
        }
        currentBlock_ = 0;
      }
    }
    while (test);
    return std::move(currentBlock_);
  }
};
} // end of namespace bpp.

#endif // _BLOCKSIZEMAFITERATOR_H_
