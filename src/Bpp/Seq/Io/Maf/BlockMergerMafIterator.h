// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _BLOCKMERGERMAFITERATOR_H_
#define _BLOCKMERGERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Merge blocks if some of their sequences are contiguous.
 *
 * The user specifies the focus species. Sequences that are not in this set will
 * be automatically merged and their coordinates removed.
 * The scores, if any, will be averaged for the block, weighted by the corresponding block sizes.
 * The pass value will be removed if it is different for the two blocks.
 * It is possible to define a maximum distance for the merging. Setting a distance of zero implies that the blocks
 * have to be exactly contiguous. Alternatively, the appropriate number of 'N' will be inserted in all species.
 * All species however have to be distant of the exact same amount.
 * Setting 'allowVariableDistance' to true will also merge in case of heterogeneous distance, adding arbitrary gaps where needed.
 */
class BlockMergerMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::vector<std::string> species_;
  std::unique_ptr<MafBlock> incomingBlock_;
  std::vector<std::string> ignoreChrs_; // These chromosomes will never be merged (ex: 'Un').
  unsigned int maxDist_;
  bool renameChimericChromosomes_;
  std::map<std::string, unsigned int> chimericChromosomeCounts_;
  bool allowVariableDistance_;

public:
  BlockMergerMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species,
      unsigned int maxDist = 0,
      bool renameChimericChromosomes = false,
      bool allowVariableDistance = false) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    incomingBlock_(nullptr),
    ignoreChrs_(),
    maxDist_(maxDist),
    renameChimericChromosomes_(renameChimericChromosomes),
    chimericChromosomeCounts_(),
    allowVariableDistance_(allowVariableDistance)
  {
    incomingBlock_ = iterator->nextBlock();
  }

private:
  BlockMergerMafIterator(const BlockMergerMafIterator& iterator) = delete;

  BlockMergerMafIterator& operator=(const BlockMergerMafIterator& iterator) = delete;

public:
  /**
   * brief Add a chromosome that should be ignored to the list.
   * @param chr The name of the chromosome to be ignored.
   */
  void ignoreChromosome(const std::string& chr)
  {
    ignoreChrs_.push_back(chr);
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _BLOCKMERGERMAFITERATOR_H_
