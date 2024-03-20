// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ORDERFILTERMAFITERATOR_H_
#define _ORDERFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator check that blocks are ordered according to a reference sequence.
 *
 * The occurrence of overlapping or unordered blocks result in an error message.
 * Alternatively, conflicting blocks can be discarded.
 */
class OrderFilterMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::string refSpecies_;
  std::string currentChr_;
  size_t previousBlockStart_;
  size_t previousBlockStop_;
  bool unsortedBlockDiscarded_;
  bool unsortedBlockThrowsException_;
  bool overlappingBlockDiscarded_;
  bool overlappingBlockThrowsException_;

public:
  /**
   * @brief Build a new OrderFilterMafIterator object.
   *
   * @param iterator The input iterator.
   * @param reference The species to use as a reference for coordinates.
   * @param unsortedBlockDiscarded Tell is unsorted blocks should be discarded
   * @param unsortedBlockThrowsException Tell is unsorted blocks should throw an exception
   * @param overlappingBlockDiscarded Tell is overlapping blocks should be discarded
   * @param overlappingBlockThrowsException Tell is overlapping blocks should throw an exception
   */
  OrderFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& reference,
      bool unsortedBlockDiscarded = true,
      bool unsortedBlockThrowsException = false,
      bool overlappingBlockDiscarded = true,
      bool overlappingBlockThrowsException = false) :
    AbstractFilterMafIterator(iterator),
    refSpecies_(reference),
    currentChr_(),
    previousBlockStart_(),
    previousBlockStop_(),
    unsortedBlockDiscarded_(unsortedBlockDiscarded),
    unsortedBlockThrowsException_(unsortedBlockThrowsException),
    overlappingBlockDiscarded_(overlappingBlockDiscarded),
    overlappingBlockThrowsException_(overlappingBlockThrowsException)
  {}

private:
  OrderFilterMafIterator(const OrderFilterMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    refSpecies_(iterator.refSpecies_),
    currentChr_(iterator.currentChr_),
    previousBlockStart_(iterator.previousBlockStart_),
    previousBlockStop_(iterator.previousBlockStop_),
    unsortedBlockDiscarded_(iterator.unsortedBlockDiscarded_),
    unsortedBlockThrowsException_(iterator.unsortedBlockThrowsException_),
    overlappingBlockDiscarded_(iterator.overlappingBlockDiscarded_),
    overlappingBlockThrowsException_(iterator.overlappingBlockThrowsException_)
  {}

  OrderFilterMafIterator& operator=(const OrderFilterMafIterator& iterator)
  {
    refSpecies_                      = iterator.refSpecies_;
    currentChr_                      = iterator.currentChr_;
    previousBlockStart_              = iterator.previousBlockStart_;
    previousBlockStop_               = iterator.previousBlockStop_;
    unsortedBlockDiscarded_          = iterator.unsortedBlockDiscarded_;
    unsortedBlockThrowsException_    = iterator.unsortedBlockThrowsException_;
    overlappingBlockDiscarded_       = iterator.overlappingBlockDiscarded_;
    overlappingBlockThrowsException_ = iterator.overlappingBlockThrowsException_;
    return *this;
  }

public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    bool testCont = true;
    while (testCont)
    {
      currentBlock_ = iterator_->nextBlock();
      if (currentBlock_)
        testCont = !parseBlock_(*currentBlock_);
      else
        testCont = false;
    }
    return move(currentBlock_);
  }

private:
  // Returns true if block is ordered with previous one
  bool parseBlock_(const MafBlock& block);
};
} // end of namespace bpp.

#endif//_ORDERFILTERMAFITERATOR_H_
