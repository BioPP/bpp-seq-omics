// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _CONCATENATEMAFITERATOR_H_
#define _CONCATENATEMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Concatenate blocks up to a certain size.
 *
 * Blocks are appended regardless of their coordinates, to form concatenated blocks of at least a given number of positions.
 * The scores, if any, will be averaged for the block, weighted by the corresponding block sizes.
 * The pass value will be removed if it is different for the blocks.
 * If a reference species is given, only block with identical chr tag will be concatenated.
 */
class ConcatenateMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::unique_ptr<MafBlock> incomingBlock_;
  unsigned int minimumSize_;
  std::string refSpecies_;

public:
  ConcatenateMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      unsigned int minimumSize,
      std::string refSpecies = "") :
    AbstractFilterMafIterator(iterator),
    incomingBlock_(nullptr),
    minimumSize_(minimumSize),
    refSpecies_(refSpecies)
  {
    incomingBlock_ = iterator->nextBlock();
  }

private:
  ConcatenateMafIterator(const ConcatenateMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    incomingBlock_(),
    minimumSize_(iterator.minimumSize_),
    refSpecies_(iterator.refSpecies_)
  {}

  ConcatenateMafIterator& operator=(const ConcatenateMafIterator& iterator)
  {
    incomingBlock_ = nullptr;
    minimumSize_ = iterator.minimumSize_;
    refSpecies_ = iterator.refSpecies_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_CONCATENATEMAFITERATOR_H_
