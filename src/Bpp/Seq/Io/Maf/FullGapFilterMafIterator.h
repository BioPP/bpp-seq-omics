//
// File: FullGapFilterMafIterator.h
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
