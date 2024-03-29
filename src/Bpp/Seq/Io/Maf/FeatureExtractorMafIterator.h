// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _FEATUREEXTRACTORMAFITERATOR_H_
#define _FEATUREEXTRACTORMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Extract alignments corresponding to sequence features given as a vector of SequenceFeature objects.
 *
 * The resulting blocks will contain the specified annotated regions.
 * Note that this iterator is not the opposite of FeatureFilterMafIterator,
 * as overlapping features will all be extracted. This iterator may therefore results
 * in duplication of original data.
 */
class FeatureExtractorMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::string refSpecies_;
  bool completeOnly_;
  bool ignoreStrand_;
  std::deque<std::unique_ptr<MafBlock>> blockBuffer_;
  std::map<std::string, RangeSet<size_t>> ranges_;

public:
  /**
   * @brief Build a new FeatureExtractor iterator.
   *
   * @param iterator The input iterator
   * @param refSpecies The reference species for feature coordinates
   * @param complete Tell if features should be extracted only if they can be extracted in full
   * @param features The set of features to extract
   * @param ignoreStrand If true, features will be extracted 'as is', without being reversed in case they are on the negative strand.
   */
  FeatureExtractorMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& refSpecies,
      const SequenceFeatureSet& features,
      bool complete = false,
      bool ignoreStrand = false) :
    AbstractFilterMafIterator(iterator),
    refSpecies_(refSpecies),
    completeOnly_(complete),
    ignoreStrand_(ignoreStrand),
    blockBuffer_(),
    ranges_()
  {
    // Build ranges:
    std::set<std::string> seqIds = features.getSequences();
    for (auto it : seqIds)
    {
      features.fillRangeCollectionForSequence(it, ranges_[it]);
    }
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _FEATUREEXTRACTORMAFITERATOR_H_
