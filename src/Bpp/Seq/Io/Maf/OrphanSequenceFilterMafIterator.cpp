//
// File: OrphanSequenceFilterMafIterator.cpp
// Authors: Julien Dutheil
// Created: Mon Mar 11 2013
//

/*
Copyright or © or Copr. Bio++ Development Team, (2013)

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

#include "OrphanSequenceFilterMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* OrphanSequenceFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_) {
    map<string, unsigned int> counts;
    for (size_t i = 0; i < currentBlock_->getNumberOfSequences(); ++i) {
      string species = currentBlock_->getSequence(i).getSpecies(); 
      counts[species]++;
    }
    bool test = counts.size() <= species_.size();
    if (test) {
      //We have to check that the species are the right one:
      bool loseCrit = false;
      bool strictCrit = true;
      bool duplicate = false;
      for (size_t i = 0; i < species_.size() && !duplicate; ++i) {
        map<string, unsigned int>::iterator it = counts.find(species_[i]);
        if (it != counts.end()) {
          loseCrit = true;
          if (rmDuplicates_ && it->second > 1) {
            //Duplicated sequences, block is discarded if asked to...
            duplicate = true;
          }
        } else {
          strictCrit = false;
        }
      }
      if (!duplicate) { //No duplicate found or duplicates are kept
        if (strictCrit || (!strict_ && loseCrit))
          return currentBlock_;
      }
    }
    //Otherwise there is at least one extra species, we get the next block...
    currentBlock_ = iterator_->nextBlock();
  }
  
  return currentBlock_;
}

