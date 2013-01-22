//
// File: FeatureFilterMafIterator.h
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

#ifndef _FEATUREFILTERMAFITERATOR_H_
#define _FEATUREFILTERMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>
#include <deque>

namespace bpp {

/**
 * @brief Remove from alignment all positions that fall within any feature from a list given as a SequenceFeatureSet object.
 *
 * Removed regions are outputed as a trash iterator.
 */
class FeatureFilterMafIterator:
  public AbstractFilterMafIterator,
  public MafTrashIterator
{
  private:
    std::string refSpecies_;
    std::deque<MafBlock*> blockBuffer_;
    std::deque<MafBlock*> trashBuffer_;
    bool keepTrashedBlocks_;
    std::map<std::string, MultiRange<size_t> > ranges_;

  public:
    FeatureFilterMafIterator(MafIterator* iterator, const std::string& refSpecies, const SequenceFeatureSet& features, bool keepTrashedBlocks) :
      AbstractFilterMafIterator(iterator),
      refSpecies_(refSpecies),
      blockBuffer_(),
      trashBuffer_(),
      keepTrashedBlocks_(keepTrashedBlocks),
      ranges_()
    {
      //Build ranges:
      std::set<std::string> seqIds = features.getSequences();
      for (std::set<std::string>::iterator it = seqIds.begin();
          it != seqIds.end();
          ++it) {
        {
          features.fillRangeCollectionForSequence(*it, ranges_[*it]);
        }
      }
    }

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

#endif //_FEATUREFILTERMAFITERATOR_H_
