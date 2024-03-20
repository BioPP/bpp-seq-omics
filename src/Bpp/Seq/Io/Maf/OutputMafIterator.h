// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _OUTPUTMAFITERATOR_H_
#define _OUTPUTMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief This iterator forward the iterator given as input after having printed its content to a file.
 */
class OutputMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  bool mask_;

public:
  OutputMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      bool mask = true) :
    AbstractFilterMafIterator(iterator),
    output_(out),
    mask_(mask)
  {
    if (output_)
      writeHeader(*output_);
  }

private:
  OutputMafIterator(const OutputMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    mask_(iterator.mask_)
  {}

  OutputMafIterator& operator=(const OutputMafIterator& iterator)
  {
    output_ = iterator.output_;
    mask_   = iterator.mask_;
    return *this;
  }

public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (output_ && currentBlock_)
      writeBlock(*output_, *currentBlock_);
    return move(currentBlock_);
  }

private:
  void writeHeader(std::ostream& out) const;
  void writeBlock(std::ostream& out, const MafBlock& block) const;
};
} // end of namespace bpp.

#endif // _OUTPUTMAFITERATOR_H_
