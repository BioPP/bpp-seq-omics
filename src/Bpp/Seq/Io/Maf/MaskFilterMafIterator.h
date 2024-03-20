// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MASKFILTERMAFITERATOR_H_
#define _MASKFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to remove regions with masked positions.
 *
 * Regions with a too high proportion of masked position in a set of species will be removed,
 * and blocks adjusted accordingly.
 */
class MaskFilterMafIterator :
  public AbstractFilterMafIterator,
  public virtual MafTrashIteratorInterface
{
private:
  std::vector<std::string> species_;
  unsigned int windowSize_;
  unsigned int step_;
  unsigned int maxMasked_;
  std::deque<std::unique_ptr<MafBlock>> blockBuffer_;
  std::deque<std::unique_ptr<MafBlock>> trashBuffer_;
  std::deque<std::vector<bool>> window_;
  bool keepTrashedBlocks_;

public:
  MaskFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species,
      unsigned int windowSize,
      unsigned int step,
      unsigned int maxMasked,
      bool keepTrashedBlocks) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    windowSize_(windowSize),
    step_(step),
    maxMasked_(maxMasked),
    blockBuffer_(),
    trashBuffer_(),
    window_(species.size()),
    keepTrashedBlocks_(keepTrashedBlocks)
  {}

public:
  std::unique_ptr<MafBlock> nextRemovedBlock()
  {
    if (trashBuffer_.size() == 0) return nullptr;
    auto block = move(trashBuffer_.front());
    trashBuffer_.pop_front();
    return block;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _MASKFILTERMAFITERATOR_H_
