//
// File: ConcatenateMafIterator.h
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

#ifndef _CONCATENATEMAFITERATOR_H_
#define _CONCATENATEMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief Concatenate blocks up to a certain size.
 *
 * Blocks are appended regardless of their coordinates, to form concatenated blocks of at least a given number of positions.
 * The scores, if any, will be averaged for the block, weighted by the corresponding block sizes.
 * The pass value will be removed if it is different for the blocks.
 */
class ConcatenateMafIterator:
  public AbstractFilterMafIterator
{
  private:
    MafBlock* incomingBlock_;
    unsigned int minimumSize_;

  public:
    ConcatenateMafIterator(MafIterator* iterator, unsigned int minimumSize) :
      AbstractFilterMafIterator(iterator),
      incomingBlock_(0),
      minimumSize_(minimumSize)
    {
      incomingBlock_ = iterator->nextBlock();
    }

  private:
    ConcatenateMafIterator(const ConcatenateMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      incomingBlock_(iterator.incomingBlock_),
      minimumSize_(iterator.minimumSize_)
    {}
    
    ConcatenateMafIterator& operator=(const ConcatenateMafIterator& iterator)
    {
      incomingBlock_ = iterator.incomingBlock_;
      minimumSize_ = iterator.minimumSize_;
      return *this;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);

};

} // end of namespace bpp.

#endif //_CONCATENATEMAFITERATOR_H_
