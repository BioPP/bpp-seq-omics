// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _FEATUREFILTERMAFITERATOR_H_
#define _FEATUREFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>
#include <memory>

namespace bpp
{
/**
 * @brief Remove from alignment all positions that fall within any feature from a list given as a SequenceFeatureSet object.
 *
 * Removed regions are outputed as a trash iterator.
 */
class FeatureFilterMafIterator :
  public AbstractFilterMafIterator,
  public MafTrashIteratorInterface
{
private:
  std::string refSpecies_;
  std::deque<std::unique_ptr<MafBlock>> blockBuffer_;
  std::deque<std::unique_ptr<MafBlock>> trashBuffer_;
  bool keepTrashedBlocks_;
  std::map<std::string, MultiRange<size_t> > ranges_;

public:
  FeatureFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& refSpecies,
      const SequenceFeatureSet& features,
      bool keepTrashedBlocks) :
    AbstractFilterMafIterator(iterator),
    refSpecies_(refSpecies),
    blockBuffer_(),
    trashBuffer_(),
    keepTrashedBlocks_(keepTrashedBlocks),
    ranges_()
  {
    // Build ranges:
    std::set<std::string> seqIds = features.getSequences();
    for (auto& it : seqIds)
    {
      {
        features.fillRangeCollectionForSequence(it, ranges_[it]);
      }
    }
  }

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

#endif//_FEATUREFILTERMAFITERATOR_H_
