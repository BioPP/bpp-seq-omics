//
// File: WindowSplitMafIterator.cpp
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

#include "WindowSplitMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

const short WindowSplitMafIterator::RAGGED_LEFT = 0;
const short WindowSplitMafIterator::RAGGED_RIGHT = 1;
const short WindowSplitMafIterator::CENTER = 2;
const short WindowSplitMafIterator::ADJUST= 3;

MafBlock* WindowSplitMafIterator::analyseCurrentBlock_() throw (Exception)
{
  while (blockBuffer_.size() == 0) {
    //Build a new series of windows:
    MafBlock* block = iterator_->nextBlock();
    if (!block) return 0; //No more block.

    size_t pos = 0;
    size_t size = windowSize_;
    size_t bSize = block->getNumberOfSites();

    switch (align_) {
      case (RAGGED_RIGHT) : { pos = bSize % windowSize_; break; }
      case (CENTER)       : { pos = (bSize % windowSize_) / 2; break; }
      case (ADJUST)       : {
          size_t x = bSize / windowSize_;
          if (x == 0) size = bSize;
          else        size = bSize / x;
          break;
        }               
      default             : { }
    }
    //cout << "Effective size: " << size << endl;
    for(size_t i = pos; i + size < bSize; i += size) {
      MafBlock* newBlock = new MafBlock();
      newBlock->setScore(block->getScore());
      newBlock->setPass(block->getPass());
      if (align_ == ADJUST) {
        if (bSize - (i + size) > 0 && bSize - (i + size) < size) {
          //cout << "Old size: " << size;
          size = bSize - i; //Adjust for last block because of rounding.
                            //this should not increase size by more than 1!
          //cout << " => new size: " << size << endl;
        }
      }
      for (size_t j = 0; j < block->getNumberOfSequences(); ++j) {
        auto_ptr<MafSequence> subseq(block->getSequence(j).subSequence(i, size));
        newBlock->addSequence(*subseq);
      }
      blockBuffer_.push_back(newBlock);
    }
    
    if (align_ == ADJUST && keepSmallBlocks_ && bSize < windowSize_) {
      blockBuffer_.push_back(block);
    } else {
      delete block;
    }
  }
 
  MafBlock* nxtBlock = blockBuffer_.front();
  blockBuffer_.pop_front();
  return nxtBlock;
}

