//
// File: DuplicateFilterMafIterator.h
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
