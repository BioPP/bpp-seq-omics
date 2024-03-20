// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _FULLGAPFILTERMAFITERATOR_H_
#define _FULLGAPFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>
#include <memory>

namespace bpp
{
/**
 * @brief Filter maf blocks to remove in each block the positions made only of gaps.
 *
 * The subset of species that should be examined is given as input. The coordinates of these
 * species will not be altered as only gap positions are removed. Other species however may be
 * altered as they might not have gap in the removed position. The coordinates for these species
 * will therefore be removed as they do not make sense anymore.
 */
class FullGapFilterMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::vector<std::string> species_;

public:
  FullGapFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species) :
    AbstractFilterMafIterator(iterator),
    species_(species)
  {}

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_FULLGAPFILTERMAFITERATOR_H_
