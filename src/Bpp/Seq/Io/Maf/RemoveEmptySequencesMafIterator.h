// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _REMOVEEMPTYSEQUENCESMAFITERATOR_H_
#define _REMOVEEMPTYSEQUENCESMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>

namespace bpp
{
/**
 * @brief Remove ga-only or unresolved and gap-only sequences in each block.
 */
class RemoveEmptySequencesMafIterator :
  public AbstractFilterMafIterator
{
private:
  bool unresolvedAsGaps_;

public:
  /**
   * @brief Creates a new RemoveEmptySequencesMafIterator object.
   *
   * @param iterator The input iterator.
   * @param unresolvedAsGaps Tell if unresolved characters (e.g. 'N') should be treated as gaps.
   */
  RemoveEmptySequencesMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      bool unresolvedAsGaps = false) :
    AbstractFilterMafIterator(iterator),
       	unresolvedAsGaps_(unresolvedAsGaps)
  {}

private:
  RemoveEmptySequencesMafIterator(const RemoveEmptySequencesMafIterator& iterator) :
    AbstractFilterMafIterator(nullptr),
    unresolvedAsGaps_(iterator.unresolvedAsGaps_)
  {}

  RemoveEmptySequencesMafIterator& operator=(const RemoveEmptySequencesMafIterator& iterator)
  {
    unresolvedAsGaps_ = iterator.unresolvedAsGaps_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_REMOVEEMPTYSEQUENCESMAFITERATOR_H_
