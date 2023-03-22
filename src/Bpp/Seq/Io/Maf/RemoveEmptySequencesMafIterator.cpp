//
// File: RemoveEmptySequencesMafIterator.cpp
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

#include "RemoveEmptySequencesMafIterator.h"

using namespace bpp;
using namespace std;

unique_ptr<MafBlock> RemoveEmptySequencesMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_)
  {
    for (size_t i = currentBlock_->getNumberOfSequences(); i > 0; --i)
    {
      const MafSequence& seq = currentBlock_->sequence(i - 1);
      bool isEmpty = true;
      if (unresolvedAsGaps_)
      {
        for (size_t j = 0; isEmpty && j < currentBlock_->getNumberOfSites(); ++j)
        {
          if (!AlphabetTools::DNA_ALPHABET->isUnresolved(seq[j]) && !AlphabetTools::DNA_ALPHABET->isGap(seq[j]))
            isEmpty = false;
        }
      }
      else
      {
        for (size_t j = 0; isEmpty && j < currentBlock_->getNumberOfSites(); ++j)
        {
          if (!AlphabetTools::DNA_ALPHABET->isGap(seq[j]))
            isEmpty = false;
        }
      }
      if (isEmpty)
      {
        currentBlock_->removeSequence(i - 1);
      }
    }
  }
  return move(currentBlock_);
}
