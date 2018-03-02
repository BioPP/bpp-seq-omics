//
// File: DuplicateFilterMafIterator.cpp
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

#include "DuplicateFilterMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* DuplicateFilterMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  while (currentBlock_) {
    bool foundRef = false;
    string chr = "";
    char strand = '+';
    size_t start = 0;
    size_t stop  = 0;
    for (size_t i = 0; i < currentBlock_->getNumberOfSequences() && !foundRef; ++i) {
      string species = currentBlock_->getSequence(i).getSpecies(); 
      if (species == ref_) {
        foundRef = true;
        chr    = currentBlock_->getSequence(i).getChromosome();
        strand = currentBlock_->getSequence(i).getStrand();
        start  = currentBlock_->getSequence(i).start();
        stop   = currentBlock_->getSequence(i).stop();
      }
    }
    if (!foundRef) {
      if (logstream_) {
        (*logstream_ << "DUPLICATE FILTER: block does not contain reference species and was removed.").endLine();
      }
      delete currentBlock_;
    } else {
      size_t occurrence = blocks_[chr][strand][start][stop]++;
      if (occurrence > 0) {
        if (logstream_) {
          (*logstream_ << "DUPLICATE FILTER: sequence in reference species was found in a previous block. New block was removed.").endLine();
        }
        delete currentBlock_;
      } else {
        return currentBlock_;
      }
    }

    //Look for the next block:
    currentBlock_ = iterator_->nextBlock();
  }
  
  return currentBlock_;
}

