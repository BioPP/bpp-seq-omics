//
// File: TableOutputMafIterator.h
// Authors: Julien Dutheil
// Created: Sun Jun 25 2017
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2017)

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

#ifndef _TABLEOUTPUTMAFITERATOR_H_
#define _TABLEOUTPUTMAFITERATOR_H_

#include "MafIterator.h"

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
  std::ostream* output_;
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
  TableOutputMafIterator(MafIterator* iterator,
                         std::ostream* out,
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
  MafBlock* analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (output_ && currentBlock_)
      writeBlock_(*output_, *currentBlock_);
    return currentBlock_;
  }

private:
  void writeBlock_(std::ostream& out, const MafBlock& block);
};
} // end of namespace bpp.

#endif//_TABLEOUTPUTMAFITERATOR_H_
