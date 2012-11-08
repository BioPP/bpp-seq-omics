//
// File: QualityFilterMafIterator.h
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

#ifndef _QUALITYFILTERMAFITERATOR_H_
#define _QUALITYFILTERMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief Filter maf blocks to remove regions with low quality.
 *
 * Regions with a too low average quality in a set of species will be removed,
 * and blocks adjusted accordingly. 
 */
class QualityFilterMafIterator:
  public AbstractFilterMafIterator,
  public MafTrashIterator
{
  private:
    std::vector<std::string> species_;
    unsigned int windowSize_;
    unsigned int step_;
    unsigned int minQual_;
    std::deque<MafBlock*> blockBuffer_;
    std::deque<MafBlock*> trashBuffer_;
    std::deque< std::vector<int> > window_;
    bool keepTrashedBlocks_;

  public:
    QualityFilterMafIterator(MafIterator* iterator, const std::vector<std::string>& species, unsigned int windowSize, unsigned int step, unsigned int minQual, bool keepTrashedBlocks) :
      AbstractFilterMafIterator(iterator),
      species_(species),
      windowSize_(windowSize),
      step_(step),
      minQual_(minQual),
      blockBuffer_(),
      trashBuffer_(),
      window_(species.size()),
      keepTrashedBlocks_(keepTrashedBlocks)
    {}

  public:
    MafBlock* nextRemovedBlock() throw (Exception) {
      if (trashBuffer_.size() == 0) return 0;
      MafBlock* block = trashBuffer_.front();
      trashBuffer_.pop_front();
      return block;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);
};

} // end of namespace bpp.

#endif //_QUALITYFILTERMAFITERATOR_H_
