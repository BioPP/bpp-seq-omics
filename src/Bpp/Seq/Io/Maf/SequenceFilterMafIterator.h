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
public:
  inline static short DUPLICATE_POLICY_DISCARD = 0;
  inline static short DUPLICATE_POLICY_SAMPLE = 1;
  inline static short DUPLICATE_POLICY_REMOVE = 2;

private:
  std::vector<std::string> species_;
  bool strict_;
  bool keep_;
  bool rmDuplicates_;
  short duplicatePolicy_;

public:
  /**
   * @param iterator The input iterator.
   * @param species The list of species names to be retained.
   * @param strict If true, then block that do not contain all species will be discarded.
   * @param keep If true, sequences not in the selection will be kept.
   * @param rmDuplicates If true, then extra filtering is dones according to the duplicate policy.
   * @param duplicatePolicy One of 
   *   DUPLICATE_POLICY_DISCARD (block that contain more than one instance for at least one species will be discarded [default for backward compatibility]),
   *   DUPLICATE_POLICY_SAMPLE (keep only one of the duplicated sequence in the block - the last one),
   *   DUPLICATE_POLICY_REMOVE (remove all duplicated sequences from the block).
   */
  SequenceFilterMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::vector<std::string>& species,
      bool strict = false,
      bool keep = false,
      bool rmDuplicates = false,
      short duplicatePolicy = DUPLICATE_POLICY_DISCARD) :
    AbstractFilterMafIterator(iterator),
    species_(species),
    strict_(strict),
    keep_(keep),
    rmDuplicates_(rmDuplicates),
    duplicatePolicy_(duplicatePolicy)
  {}

private:
  SequenceFilterMafIterator(const SequenceFilterMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    species_(iterator.species_),
    strict_(iterator.strict_),
    keep_(iterator.keep_),
    rmDuplicates_(iterator.rmDuplicates_),
    duplicatePolicy_(iterator.duplicatePolicy_)
  {}

  SequenceFilterMafIterator& operator=(const SequenceFilterMafIterator& iterator)
  {
    species_         = iterator.species_;
    strict_          = iterator.strict_;
    keep_            = iterator.keep_;
    rmDuplicates_    = iterator.rmDuplicates_;
    duplicatePolicy_ = iterator.duplicatePolicy_;
    return *this;
  }

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif // _SEQUENCEFILTERMAFITERATOR_H_
