// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _DUPLICATEFILTERMAFITERATOR_H_
#define _DUPLICATEFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to remove duplicated blocks, according to a reference sequence).
 */
class DuplicateFilterMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::string ref_;
  /**
   * Contains the list of 'seen' block, as [chr][strand][start][stop]
   */
  std::map< std::string, std::map< char, std::map< size_t, std::map< size_t, size_t > > > > blocks_;

public:
  /**
   * @param iterator The input iterator.
   * @param reference The reference species name.
   */
  DuplicateFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& reference) :
    AbstractFilterMafIterator(iterator),
    ref_(reference),
    blocks_()
  {}

private:
  DuplicateFilterMafIterator(const DuplicateFilterMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    ref_(iterator.ref_),
    blocks_(iterator.blocks_)
  {}

  DuplicateFilterMafIterator& operator=(const DuplicateFilterMafIterator& iterator)
  {
    ref_    = iterator.ref_;
    blocks_ = iterator.blocks_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_DUPLICATEFILTERMAFITERATOR_H_
