//
// File: OutputAlignmentMafIterator.cpp
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

#include "OutputAlignmentMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* OutputAlignmentMafIterator::analyseCurrentBlock_() throw (Exception)
{
  MafBlock* block = iterator_->nextBlock();
  if (block) {
    if (output_) {
      writeBlock(*output_, *block);
    } else {
      string file = file_;
      TextTools::replaceAll(file, "%i", TextTools::toString(++currentBlockIndex_));
      std::ofstream output(file.c_str(), ios::out);
      writeBlock(output, *block);
    }
  }
  return block;
}

void OutputAlignmentMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const {
  //First get alignment:
  AlignedSequenceContainer aln(&AlphabetTools::DNA_ALPHABET);
  //We cannot copy directly the container because we want to convert from MafSequence to BasicSequence (needed for renaiming):
  SequenceContainerTools::convertContainer<AlignedSequenceContainer, AlignedSequenceContainer, BasicSequence>(block.getAlignment(), aln);
  //Format sequence names:
  vector<string> names(aln.getNumberOfSequences());
  for (unsigned int i = 0; i < aln.getNumberOfSequences(); ++i) {
    const MafSequence& mafseq = block.getSequence(i);
    names[i] = mafseq.getSpecies() + "-" + mafseq.getChromosome() + "(" + mafseq.getStrand() + ")/" + TextTools::toString(mafseq.start() + 1) + "-" + TextTools::toString(mafseq.stop() + 1);
  }
  aln.setSequencesNames(names);
  writer_->writeAlignment(out, aln);
}

