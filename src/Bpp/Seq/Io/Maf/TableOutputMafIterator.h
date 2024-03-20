// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _TABLEOUTPUTMAFITERATOR_H_
#define _TABLEOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator outputs sequence states for selected species and positions
 */
class TableOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::vector<std::string> species_;
  std::string refSpecies_;

public:
  /**
   * @brief Build a new TableOutputMafIterator object.
   *
   * @param iterator The input iterator.
   * @param out The output stream where to write the table file.
   * @param species A list of species for which sequence content should be output (one column per selected species).
   * In case one species is duplicated in a block, the first sequence will be used.
   * @param reference The species to use as a reference for coordinates.
   * It does not have to be one of the selected species for output.
   */
  TableOutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::vector<std::string>& species,
      const std::string& reference) :
    AbstractFilterMafIterator(iterator),
    output_(out), species_(species), refSpecies_(reference)
  {
    // Write header:
    *output_ << "Chromosome\tPosition";
    for (const std::string& sp : species)
    {
      *output_ << "\t" << sp;
    }
    *output_ << std::endl;
  }

private:
  TableOutputMafIterator(const TableOutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    species_(iterator.species_),
    refSpecies_(iterator.refSpecies_)
  {}

  TableOutputMafIterator& operator=(const TableOutputMafIterator& iterator)
  {
    output_ = iterator.output_;
    species_ = iterator.species_;
    refSpecies_ = iterator.refSpecies_;
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

#endif//_TABLEOUTPUTMAFITERATOR_H_
