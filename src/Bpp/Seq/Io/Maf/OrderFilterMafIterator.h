//
// File: OrderFilterMafIterator.h
// Authors: Julien Dutheil
// Created: Sat Jan 20 2018
//

/*
Copyright or © or Copr. Bio++ Development Team, (2018)

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

#ifndef _ORDERFILTERMAFITERATOR_H_
#define _ORDERFILTERMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief This iterator check that blocks are ordered according to a reference sequence.
 *
 * The occurrence of overlapping or unordered blocks result in an error message.
 * Alternatively, conflicting blocks can be discarded.
 */
class OrderFilterMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::string refSpecies_;
    std::string currentChr_;
    size_t previousBlockStart_;
    size_t previousBlockStop_;
    bool unsortedBlockDiscarded_;
    bool unsortedBlockThrowsException_;
    bool overlappingBlockDiscarded_;
    bool overlappingBlockThrowsException_;

  public:
    /**
     * @brief Build a new OrderFilterMafIterator object.
     *
     * @param iterator The input iterator.
     * @param reference The species to use as a reference for coordinates.
     * @param unsortedBlockDiscarded Tell is unsorted blocks should be discarded
     * @param unsortedBlockThrowsException Tell is unsorted blocks should throw an exception
     * @param overlappingBlockDiscarded Tell is overlapping blocks should be discarded
     * @param overlappingBlockThrowsException Tell is overlapping blocks should throw an exception
     */
    OrderFilterMafIterator(MafIterator* iterator,
        const std::string& reference,
        bool unsortedBlockDiscarded = true,
        bool unsortedBlockThrowsException = false,
        bool overlappingBlockDiscarded = true,
        bool overlappingBlockThrowsException = false) :
      AbstractFilterMafIterator(iterator),
      refSpecies_(reference),
      currentChr_(),
      previousBlockStart_(),
      previousBlockStop_(),
      unsortedBlockDiscarded_(unsortedBlockDiscarded),
      unsortedBlockThrowsException_(unsortedBlockThrowsException),
      overlappingBlockDiscarded_(overlappingBlockDiscarded),
      overlappingBlockThrowsException_(overlappingBlockThrowsException)
    {}

  private:
    OrderFilterMafIterator(const OrderFilterMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      refSpecies_(iterator.refSpecies_),
      currentChr_(iterator.currentChr_),
      previousBlockStart_(iterator.previousBlockStart_),
      previousBlockStop_(iterator.previousBlockStop_),
      unsortedBlockDiscarded_(iterator.unsortedBlockDiscarded_),
      unsortedBlockThrowsException_(iterator.unsortedBlockThrowsException_),
      overlappingBlockDiscarded_(iterator.overlappingBlockDiscarded_),
      overlappingBlockThrowsException_(iterator.overlappingBlockThrowsException_)
    {}
    
    OrderFilterMafIterator& operator=(const OrderFilterMafIterator& iterator)
    {
      refSpecies_                      = iterator.refSpecies_;
      currentChr_                      = iterator.currentChr_;
      previousBlockStart_              = iterator.previousBlockStart_;
      previousBlockStop_               = iterator.previousBlockStop_;
      unsortedBlockDiscarded_          = iterator.unsortedBlockDiscarded_;
      unsortedBlockThrowsException_    = iterator.unsortedBlockThrowsException_;
      overlappingBlockDiscarded_       = iterator.overlappingBlockDiscarded_;
      overlappingBlockThrowsException_ = iterator.overlappingBlockThrowsException_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_() {
      bool testCont = true;
      while (testCont) {
        currentBlock_ = iterator_->nextBlock();
        if (currentBlock_)
          testCont = !parseBlock_(*currentBlock_);
        else
          testCont = false;
      }
      return currentBlock_;
    }

  private:
    //Returns true if block is ordered with previous one
    bool parseBlock_(const MafBlock& block);
};

} // end of namespace bpp.

#endif //_ORDERFILTERMAFITERATOR_H_
