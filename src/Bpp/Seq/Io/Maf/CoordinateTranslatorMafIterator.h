// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _COORDINATETRANSLATORMAFITERATOR_H_
#define _COORDINATETRANSLATORMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

// From bpp-core:
#include <Bpp/Numeric/DataTable.h>

namespace bpp
{
/**
 * @brief Translate features coordinates from one species to another, based on the alignment.
 *
 * This filter is similar in principle to the UCSC "liftOver" utility and software alike.
 * For now, only write a text file with all coordinates from reference and corresponding target sequence.
 */
class CoordinateTranslatorMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::string referenceSpecies_;
  std::string targetSpecies_;
  std::map<std::string, SequenceFeatureSet*> inputFeaturesPerChr_;
  std::ostream& output_;
  bool outputClosestCoordinate_;

public:
  /**
   * @brief Build a new CoordinateTranslator iterator.
   *
   * @param iterator The input iterator
   * @param referenceSpecies The reference species for feature coordinates
   * @param targetSpecies The target species for which features coordinates should be translated
   * @param features The set of features to lift over
   * @param output Output stream for translated coordinates
   * @param outputClosestCoordinate In case the target sequence has a gap at the corresponding position,
   *        tells if the previous non-gap position should be returned, or NA.
   */
  CoordinateTranslatorMafIterator(
    std::shared_ptr<MafIteratorInterface> iterator,
    const std::string& referenceSpecies,
    const std::string& targetSpecies,
    const SequenceFeatureSet& features,
    std::ostream& output,
    bool outputClosestCoordinate = true) :
    AbstractFilterMafIterator(iterator),
    referenceSpecies_(referenceSpecies),
    targetSpecies_(targetSpecies),
    inputFeaturesPerChr_(),
    output_(output),
    outputClosestCoordinate_(outputClosestCoordinate)
  {
    // Sort features per chromosome for a faster access:
    auto seqIds = features.getSequences();
    for (auto it : seqIds)
    {
      {
        inputFeaturesPerChr_[it] = features.getSubsetForSequence(it);
      }
    }
    output_ << "chr.ref\tstrand.ref\tbegin.ref\tend.ref\tchr.target\tstrand.target\tbegin.target\tend.target" << std::endl;
  }

  virtual ~CoordinateTranslatorMafIterator()
  {
    // Clean sorted features.
    for (std::map<std::string, SequenceFeatureSet*>::iterator it = inputFeaturesPerChr_.begin();
         it != inputFeaturesPerChr_.end();
         it++)
    {
      delete it->second;
    }
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_COORDINATETRANSLATORMAFITERATOR_H_
