//
// File: FeatureExtractorMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Sep 07 2010
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2010)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _FEATUREEXTRACTORMAFITERATOR_H_
#define _FEATUREEXTRACTORMAFITERATOR_H_

#include "MafIterator.h"

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
  std::deque<MafBlock*> blockBuffer_;
  std::map<std::string, RangeSet<size_t> > ranges_;

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
  FeatureExtractorMafIterator(MafIterator* iterator, const std::string& refSpecies, const SequenceFeatureSet& features, bool complete = false, bool ignoreStrand = false) :
    AbstractFilterMafIterator(iterator),
    refSpecies_(refSpecies),
    completeOnly_(complete),
    ignoreStrand_(ignoreStrand),
    blockBuffer_(),
    ranges_()
  {
    // Build ranges:
    std::set<std::string> seqIds = features.getSequences();
    for (std::set<std::string>::iterator it = seqIds.begin();
         it != seqIds.end();
         ++it)
    {
      {
        features.fillRangeCollectionForSequence(*it, ranges_[*it]);
      }
    }
  }

private:
  MafBlock* analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_FEATUREEXTRACTORMAFITERATOR_H_
