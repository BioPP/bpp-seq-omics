//
// File: CoordinatesOutputMafIterator.cpp
// Authors: Julien Dutheil
// Created: Mon Jun 02 2014
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

#include "CoordinatesOutputMafIterator.h"

using namespace bpp;
using namespace std;

void CoordinatesOutputMafIterator::writeHeader_(ostream& out) const
{
  for (size_t i = 0; i < species_.size(); ++i) {
    if (i > 0) out << "\t";
    string sp = species_[i];
    out << sp << ".chr\t" << sp << ".strand\t" << sp << ".start\t" << sp << ".stop";
    if (includeSrcSize_)
      out << "\t" << sp << ".src";
  }
  out << endl;
}

MafBlock* CoordinatesOutputMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_) {
    for (size_t i = 0; i < species_.size(); ++i) {
      if (i > 0) *output_ << "\t";
      vector<const MafSequence*> seqs = currentBlock_->getSequencesForSpecies(species_[i]);
      if (seqs.size() > 1)
        throw Exception("CoordinatesOutputMafIterator::analyseCurrentBlock_(). There is more than one sequence for species '" + species_[i] + "' in current block.");
      else if (seqs.size() == 0) {
        *output_ << "NA\tNA\tNA\tNA";
        if (includeSrcSize_)
          *output_ << "\tNA";
      } else {
        *output_ << seqs[0]->getChromosome() << "\t" << seqs[0]->getStrand() << "\t" << seqs[0]->start() << "\t" << seqs[0]->stop();
        if (includeSrcSize_)
          *output_ << "\t" << seqs[0]->getSrcSize();
      }
    }
    *output_ << endl;
  }
  return currentBlock_;
}

