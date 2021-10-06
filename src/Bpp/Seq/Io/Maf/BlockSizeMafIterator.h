//
// File: BlockSizeMafIterator.h
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

#ifndef _BLOCKSIZEMAFITERATOR_H_
#define _BLOCKSIZEMAFITERATOR_H_

#include "MafIterator.h"

// From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp
{
/**
 * @brief Filter maf blocks to keep only the ones with a minimum number of species.
 */
class BlockSizeMafIterator :
  public AbstractFilterMafIterator
{
private:
  unsigned int minSize_;

public:
  BlockSizeMafIterator(MafIterator* iterator, unsigned int minSize) :
    AbstractFilterMafIterator(iterator),
    minSize_(minSize)
  {}

private:
  MafBlock* analyseCurrentBlock_() throw (Exception)
  {
    bool test;
    do
    {
      currentBlock_ = iterator_->nextBlock();
      if (!currentBlock_) break;
      test = (currentBlock_->getNumberOfSequences() < minSize_);
      if (test)
      {
        if (logstream_)
        {
          (*logstream_ << "BLOCK SIZE FILTER: block " << currentBlock_->getDescription() << " with size " << currentBlock_->getNumberOfSites() << " was discarded.").endLine();
        }
        delete currentBlock_;
        currentBlock_ = 0;
      }
    }
    while (test);
    return currentBlock_;
  }
};
} // end of namespace bpp.

#endif//_BLOCKSIZEMAFITERATOR_H_
