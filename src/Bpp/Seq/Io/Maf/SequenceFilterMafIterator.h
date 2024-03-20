// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _SEQUENCEFILTERMAFITERATOR_H_
#define _SEQUENCEFILTERMAFITERATOR_H_

#include "AbstractMafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to keep a subset of sequences, given their name.
 * this filter can also be used to filter block which contain a certain set of sequences.
 *
 * Typical usage:
 * - strict=yes, keep=no: extract the species from the list, only if all of them are present in a block.
 * - strict=no, keep=no: extract the species from the list, at least the one which are there.
 * - strict=yes, keep=yes: filter blocks to retain only the ones that contain at least all species from the list.
 * Blocks that are empty after the filtering are removed.
 */
class SequenceFilterMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::vector<std::string> species_;
  bool strict_;
  bool keep_;
  bool rmDuplicates_;

public:
  /**
   * @param iterator The input iterator.
   * @param species The list of species names to be retained.
   * @param strict If true, then block that do not contain all species will be discarded.
   * @param keep If true, sequences not in the selection will be kept.
   * @param rmDuplicates If true, block that contain more than one instance for at least one species will be discarded.
   */
  SequenceFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species,
      bool strict = false,
      bool keep = false,
      bool rmDuplicates = false) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    strict_(strict),
    keep_(keep),
    rmDuplicates_(rmDuplicates)
  {}

private:
  SequenceFilterMafIterator(const SequenceFilterMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    species_(iterator.species_),
    strict_(iterator.strict_),
    keep_(iterator.keep_),
    rmDuplicates_(iterator.rmDuplicates_)
  {}

  SequenceFilterMafIterator& operator=(const SequenceFilterMafIterator& iterator)
  {
    species_       = iterator.species_;
    strict_        = iterator.strict_;
    keep_          = iterator.keep_;
    rmDuplicates_  = iterator.rmDuplicates_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_SEQUENCEFILTERMAFITERATOR_H_
