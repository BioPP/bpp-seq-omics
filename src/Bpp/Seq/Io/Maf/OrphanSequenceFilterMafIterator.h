// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ORPHANSEQUENCEFILTERMAFITERATOR_H_
#define _ORPHANSEQUENCEFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to keep a the ones which display a specified combination of species.
 *
 * This filter is typically used to retrieve "orphan" sequences, that is sequences only present in one (set of) species.
 */
class OrphanSequenceFilterMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::vector<std::string> species_;
  bool strict_;
  bool rmDuplicates_;

public:
  /**
   * @param iterator The input iterator.
   * @param species The list of species names to be retained.
   * @param strict If true, then block that do not contain all of the specified species will be discarded.
   * @param keep If true, then additional species sequences will be kept.
   * @param rmDuplicates If true, block that contain more than one instance for at least one species will be discarded.
   */
  OrphanSequenceFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species,
      bool strict = false,
      bool keep = false,
      bool rmDuplicates = false) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    strict_(strict),
    rmDuplicates_(rmDuplicates)
  {}

private:
  OrphanSequenceFilterMafIterator(const OrphanSequenceFilterMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    species_(iterator.species_),
    strict_(iterator.strict_),
    rmDuplicates_(iterator.rmDuplicates_)
  {}

  OrphanSequenceFilterMafIterator& operator=(const OrphanSequenceFilterMafIterator& iterator)
  {
    species_       = iterator.species_;
    strict_        = iterator.strict_;
    rmDuplicates_  = iterator.rmDuplicates_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _ORPHANSEQUENCEFILTERMAFITERATOR_H_
