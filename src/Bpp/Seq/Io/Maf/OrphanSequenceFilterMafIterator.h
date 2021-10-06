//
// File: OrphanSequenceFilterMafIterator.h
// Authors: Julien Dutheil
// Created: Mon Mar 11 2013
//

/*
   Copyright or © or Copr. Bio++ Development Team, (2013)

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

#ifndef _ORPHANSEQUENCEFILTERMAFITERATOR_H_
#define _ORPHANSEQUENCEFILTERMAFITERATOR_H_

#include "MafIterator.h"

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
  OrphanSequenceFilterMafIterator(MafIterator* iterator,
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
  MafBlock* analyseCurrentBlock_();
};
} // end of namespace bpp.

#endif//_ORPHANSEQUENCEFILTERMAFITERATOR_H_
