//
// File: RemoveEmptySequencesMafIterator.h
// Authors: Julien Dutheil
// Created: Tue Apr 26 2016
//

/*
Copyright or © or Copr. Bio++ Development Team, (2014)

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

#ifndef _REMOVEEMPTYSEQUENCESMAFITERATOR_H_
#define _REMOVEEMPTYSEQUENCESMAFITERATOR_H_

#include "MafIterator.h"

//From the STL:
#include <iostream>
#include <string>

namespace bpp {

/**
 * @brief Remove ga-only or unresolved and gap-only sequences in each block.
 */
class RemoveEmptySequencesMafIterator:
  public AbstractFilterMafIterator
{
  private:
    bool unresolvedAsGaps_;

  public:
    /**
     * @brief Creates a new RemoveEmptySequencesMafIterator object.
     *
     * @param iterator The input iterator.
     * @param unresolvedAsGaps Tell if unresolved characters (e.g. 'N') should be treated as gaps.
     */
    RemoveEmptySequencesMafIterator(MafIterator* iterator, bool unresolvedAsGaps = false):
      AbstractFilterMafIterator(iterator), unresolvedAsGaps_(unresolvedAsGaps)
    {}

  private:
    RemoveEmptySequencesMafIterator(const RemoveEmptySequencesMafIterator& iterator) :
      AbstractFilterMafIterator(nullptr),
      unresolvedAsGaps_(iterator.unresolvedAsGaps_)
    {}
    
    RemoveEmptySequencesMafIterator& operator=(const RemoveEmptySequencesMafIterator& iterator)
    {
      unresolvedAsGaps_ = iterator.unresolvedAsGaps_;
      return *this;
    }

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception);

};

} // end of namespace bpp.

#endif //_REMOVEEMPTYSEQUENCESMAFITERATOR_H_
