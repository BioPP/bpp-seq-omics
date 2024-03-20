// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _SEQUENCELDHOTOUTPUTMAFITERATOR_H_
#define _SEQUENCELDHOTOUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

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
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& file,
      bool completeOnly = true,
      const std::string& reference = "") :
    AbstractFilterMafIterator(iterator),
    file_(file),
    refSpecies_(reference),
    currentBlockIndex_(0),
    completeOnly_(completeOnly)
  {}

  virtual ~SequenceLDhotOutputMafIterator() {}

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
  std::unique_ptr<MafBlock> analyseCurrentBlock_();

  void writeBlock(std::ostream& out, const MafBlock& block) const;
};
} // end of namespace bpp.

#endif // _SEQUENCEOUTPUTLDHOTMAFITERATOR_H_
