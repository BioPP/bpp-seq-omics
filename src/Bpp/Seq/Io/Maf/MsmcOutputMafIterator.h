// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MSMCOUTPUTMAFITERATOR_H_
#define _MSMCOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator outputs all SNPs in the format readable by MSMC
 *
 * See https://github.com/stschiff/msmc for a format description.
 */
class MsmcOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::vector<std::string> species_;
  std::string refSpecies_;
  std::string currentChr_;
  size_t lastPosition_;
  unsigned int nbOfCalledSites_;

public:
  /**
   * @brief Build a new MsmcOutputMafIterator object.
   *
   * @param iterator The input iterator.
   * @param out The output stream where to write the MSMC file.
   * @param species A list of at least two species to compute SNPs.
   * Only blocks containing at least these two species will be used.
   * In case one species is duplicated in a block, the first sequence will be used.
   * @param reference The species to use as a reference for coordinates.
   * It does not have to be one of the selected species on which SNPs are computed.
   */
  MsmcOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::vector<std::string>& species,
      const std::string& reference) :
    AbstractFilterMafIterator(iterator),
    output_(out), species_(species), refSpecies_(reference),
    currentChr_(""), lastPosition_(0), nbOfCalledSites_(0)
  {}

private:
  MsmcOutputMafIterator(const MsmcOutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    species_(iterator.species_),
    refSpecies_(iterator.refSpecies_),
    currentChr_(iterator.currentChr_),
    lastPosition_(iterator.lastPosition_),
    nbOfCalledSites_(iterator.nbOfCalledSites_)
  {}

  MsmcOutputMafIterator& operator=(const MsmcOutputMafIterator& iterator)
  {
    output_ = iterator.output_;
    species_ = iterator.species_;
    refSpecies_ = iterator.refSpecies_;
    currentChr_ = iterator.currentChr_;
    lastPosition_ = iterator.lastPosition_;
    nbOfCalledSites_ = iterator.nbOfCalledSites_;
    return *this;
  }

public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (output_ && currentBlock_)
      writeBlock_(*output_, *currentBlock_);
    return move(currentBlock_);
  }

private:
  void writeBlock_(std::ostream& out, const MafBlock& block);
};
} // end of namespace bpp.

#endif // _MSMCOUTPUTMAFITERATOR_H_
