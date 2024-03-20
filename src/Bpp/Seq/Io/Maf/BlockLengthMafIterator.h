// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _BLOCKLENGTHMAFITERATOR_H_
#define _BLOCKLENGTHMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to keep only the ones with a minimum number of sites.
 */
class BlockLengthMafIterator :
  public AbstractFilterMafIterator
{
private:
  size_t minLength_;

public:
  BlockLengthMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      size_t minLength) :
    AbstractFilterMafIterator(iterator),
    minLength_(minLength)
  {}

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_() override
  {
    bool test;
    do
    {
      currentBlock_ = iterator_->nextBlock();
      if (!currentBlock_) break;
      test = (currentBlock_->getNumberOfSites() < minLength_);
      if (test)
      {
        if (logstream_)
        {
          (*logstream_ << "BLOCK LENGTH FILTER: block " << currentBlock_->getDescription() << " with size " << currentBlock_->getNumberOfSites() << " was discarded.").endLine();
        }
        currentBlock_ = 0;
      }
    }
    while (test);
    return std::move(currentBlock_);
  }
};
} // end of namespace bpp.

#endif // _BLOCKLENGTHMAFITERATOR_H_
