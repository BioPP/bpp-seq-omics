// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ENTROPYFILTERMAFITERATOR_H_
#define _ENTROPYFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks highly divergent regions.
 *
 * This iterators takes two parameters: g=maxEnt and n=maxPos. Windows with more than n positions with a entropy higher than maxEnt will be discarded.
 * In addition, consecutives patterns are only counted once.
 * In case a sequence from the list is missing, it can be either ignored or counted as a full sequence of gaps.
 */
class EntropyFilterMafIterator :
  public AbstractFilterMafIterator,
  public virtual MafTrashIteratorInterface
{
private:
  std::vector<std::string> species_;
  unsigned int windowSize_;
  unsigned int step_;
  double maxEnt_;
  unsigned int maxPos_;
  std::deque<std::unique_ptr<MafBlock>> blockBuffer_;
  std::deque<std::unique_ptr<MafBlock>> trashBuffer_;
  std::deque<unsigned int> window_;
  bool keepTrashedBlocks_;
  bool missingAsGap_;
  bool ignoreGaps_;

public:
  EntropyFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species,
      unsigned int windowSize,
      unsigned int step,
      double maxEnt,
      unsigned int maxPos,
      bool keepTrashedBlocks,
      bool missingAsGap,
      bool ignoreGaps) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    windowSize_(windowSize),
    step_(step),
    maxEnt_(maxEnt),
    maxPos_(maxPos),
    blockBuffer_(),
    trashBuffer_(),
    window_(species.size()),
    keepTrashedBlocks_(keepTrashedBlocks),
    missingAsGap_(missingAsGap),
    ignoreGaps_(ignoreGaps)
  {}

public:
  std::unique_ptr<MafBlock> nextRemovedBlock()
  {
    if (trashBuffer_.size() == 0) return 0;
    auto block = std::move(trashBuffer_.front());
    trashBuffer_.pop_front();
    return block;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _ENTROPYFILTERMAFITERATOR_H_
