//
// File: SequenceFilterMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Sep 07 2010
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
