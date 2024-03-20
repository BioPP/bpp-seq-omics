// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "WindowSplitMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

const short WindowSplitMafIterator::RAGGED_LEFT = 0;
const short WindowSplitMafIterator::RAGGED_RIGHT = 1;
const short WindowSplitMafIterator::CENTER = 2;
const short WindowSplitMafIterator::ADJUST = 3;

unique_ptr<MafBlock> WindowSplitMafIterator::analyseCurrentBlock_()
{
  // Note 02/08/21: For now overlapping windows are only supported with the RAGGED_LEFT option. It could be easily generalized for cases where the window size is a multiple of the step value. More general cases are more tricky to implement, in particular for the ADJUST case.
  while (blockBuffer_.size() == 0)
  {
    // Build a new series of windows:
    auto block = iterator_->nextBlock();
    if (!block)
      return 0; // No more block.

    size_t pos = 0;
    size_t size = windowSize_;
    size_t bSize = block->getNumberOfSites();

    switch (align_)
    {
    case (RAGGED_RIGHT): { pos = bSize % windowSize_; windowStep_ = size; break; }
    case (CENTER): { pos = (bSize % windowSize_) / 2; windowStep_ = size; break; }
    case (ADJUST): {
      size_t x = bSize / windowSize_;
      if (x > 0)
        size = bSize / x;
      windowStep_ = size;
      break;
    }
    default: { }
    }
    // cout << "Effective size: " << size << endl;
    for (size_t i = pos; i + size <= bSize; i += windowStep_)
    {
      auto newBlock = make_unique<MafBlock>();
      newBlock->setScore(block->getScore());
      newBlock->setPass(block->getPass());
      if (align_ == ADJUST)
      {
        if (bSize - (i + size) > 0 && bSize - (i + size) < size)
        {
          // cout << "Old size: " << size;
          size = bSize - i; // Adjust for last block because of rounding.
                            // this should not increase size by more than 1!
          // cout << " => new size: " << size << endl;
        }
      }
      for (size_t j = 0; j < block->getNumberOfSequences(); ++j)
      {
        auto subseq = block->sequence(j).subSequence(i, size);
        newBlock->addSequence(subseq);
      }
      blockBuffer_.push_back(std::move(newBlock));
    }

    if (align_ == ADJUST && keepSmallBlocks_ && bSize < windowSize_)
    {
      blockBuffer_.push_back(std::move(block));
    }
  }

  auto nxtBlock = std::move(blockBuffer_.front());
  blockBuffer_.pop_front();
  return nxtBlock;
}
