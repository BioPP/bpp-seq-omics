//
// File: SequenceLDhotOutputMafIterator.h
// Authors: Julien Dutheil
// Created: Thr May 12 2016
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

#ifndef _SEQUENCELDHOTOUTPUTMAFITERATOR_H_
#define _SEQUENCELDHOTOUTPUTMAFITERATOR_H_

#include "MafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/Io/OSequence.h>

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator forward the iterator given as input after having printed its content to an alignment file.
 */
class SequenceLDhotOutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::string file_;
  std::string refSpecies_;
  unsigned int currentBlockIndex_;
  bool completeOnly_;

public:
  /**
   * @brief Creates a SequenceLDhotOutputMafIterator
   *
   * @param iterator The input iterator
   * @param file A string describing the path to the output files. Each block will be written to a distinct file.
   * @param completeOnly Only export complete sites (no gap, no unresolved character)
   * If "file" is a fixed string, it will only contain the last block. Using the %i code in the file name allows to generate one file per block, %i denoting the block index.
   * @param reference [optional] specify a reference species which can be used to configure file names
   * (for instance using coordinates information).
   */
  SequenceLDhotOutputMafIterator(
    MafIterator* iterator,
    const std::string& file,
    bool completeOnly = true,
    const std::string& reference = "") :
    AbstractFilterMafIterator(iterator),
    file_(file),
    refSpecies_(reference),
    currentBlockIndex_(0),
    completeOnly_(completeOnly)
  {}

  ~SequenceLDhotOutputMafIterator() {}

private:
  SequenceLDhotOutputMafIterator(const SequenceLDhotOutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    file_(iterator.file_),
    refSpecies_(iterator.refSpecies_),
    currentBlockIndex_(iterator.currentBlockIndex_),
    completeOnly_(iterator.completeOnly_)
  {}

  SequenceLDhotOutputMafIterator& operator=(const SequenceLDhotOutputMafIterator& iterator)
  {
    file_ = iterator.file_;
    refSpecies_ = iterator.refSpecies_;
    currentBlockIndex_ = iterator.currentBlockIndex_;
    completeOnly_ = iterator.completeOnly_;
    return *this;
  }

private:
  MafBlock* analyseCurrentBlock_();

  void writeBlock(std::ostream& out, const MafBlock& block) const;
};
} // end of namespace bpp.

#endif//_SEQUENCEOUTPUTLDHOTMAFITERATOR_H_
