//
// File: BlockMergerMafIterator.cpp
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

#include "BlockMergerMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* BlockMergerMafIterator::analyseCurrentBlock_() throw (Exception)
{
  if (!incomingBlock_) return 0;
  currentBlock_  = incomingBlock_;
  incomingBlock_ = iterator_->nextBlock();
  while (incomingBlock_) {
    size_t globalSpace = 0;
    for (size_t i = 0; i < species_.size(); ++i) {
      try {
        const MafSequence* seq1 = &currentBlock_->getSequenceForSpecies(species_[i]); 
        const MafSequence* seq2 = &incomingBlock_->getSequenceForSpecies(species_[i]);
        if (!seq1->hasCoordinates() || !seq2->hasCoordinates())
          throw Exception("BlockMergerMafIterator::nextBlock. Species '" + species_[i] + "' is missing coordinates in at least one block.");

        if (seq1->stop() >= seq2->start())
          return currentBlock_;
        size_t space = seq2->start() - seq1->stop() - 1;
        if (space > maxDist_)
          return currentBlock_;
        if (i == 0)
          globalSpace = space;
        else {
          if (space != globalSpace)
            return currentBlock_;
        }
        if (seq1->getChromosome() != seq2->getChromosome()
         || VectorTools::contains(ignoreChrs_, seq1->getChromosome())
         || VectorTools::contains(ignoreChrs_, seq2->getChromosome())
         || seq1->getStrand() != seq2->getStrand()
         || seq1->getSrcSize() != seq2->getSrcSize())
        {
          //There is a syntheny break in this sequence, so we do not merge the blocks.
          return currentBlock_;
        }
      } catch (SequenceNotFoundException& snfe) {
        //At least one block does not contain the sequence.
        //We don't merge the blocks:
        return currentBlock_;
      }
    }
    //We merge the two blocks:
    if (logstream_) {
      (*logstream_ << "BLOCK MERGER: merging two consecutive blocks.").endLine();
    }
    vector<string> sp1 = currentBlock_->getSpeciesList();
    vector<string> sp2 = incomingBlock_->getSpeciesList();
    vector<string> allSp = VectorTools::unique(VectorTools::vectorUnion(sp1, sp2));
    //We need to create a new MafBlock:
    MafBlock* mergedBlock = new MafBlock();
    //We average the score and pass values:
    unsigned int p1 = currentBlock_->getPass();
    unsigned int p2 = incomingBlock_->getPass();
    if (p1 == p2) mergedBlock->setPass(p1);
    double s1 = currentBlock_->getScore();
    double n1 = static_cast<double>(currentBlock_->getNumberOfSites());
    double s2 = incomingBlock_->getScore();
    double n2 = static_cast<double>(incomingBlock_->getNumberOfSites());
    mergedBlock->setScore((s1 * n1 + s2 * n2) / (n1 + n2));

    //Now fill the new block:
    for (size_t i = 0; i < allSp.size(); ++i) {
      auto_ptr<MafSequence> seq;
      try {
        seq.reset(new MafSequence(currentBlock_->getSequenceForSpecies(allSp[i])));

        //Check is there is a second sequence:
        try {
          auto_ptr<MafSequence> tmp(new MafSequence(incomingBlock_->getSequenceForSpecies(allSp[i])));
          string ref1 = seq->getDescription(), ref2 = tmp->getDescription();
          //Add spacer if needed:
          if (globalSpace > 0) {
            if (logstream_) {
              (*logstream_ << "BLOCK MERGER: a spacer of size " << globalSpace <<" is inserted in sequence for species " << allSp[i] << ".").endLine();
            }
            seq->append(vector<int>(globalSpace, AlphabetTools::DNA_ALPHABET.getUnknownCharacterCode()));
          }
          if (seq->getChromosome() != tmp->getChromosome()) {
            seq->setChromosome(seq->getChromosome() + "-" + tmp->getChromosome());
            seq->removeCoordinates();
          }
          if (seq->getStrand() != tmp->getStrand()) {
            seq->setStrand('?');
            seq->removeCoordinates();
          }
          if (seq->getName() != tmp->getName())
            tmp->setName(seq->getName()); //force name conversion to prevent exception in 'merge'.
          seq->merge(*tmp);
          if (logstream_) {
            (*logstream_ << "BLOCK MERGER: merging " << ref1 << " with " << ref2 << " into " << seq->getDescription()).endLine();
          }
        } catch (SequenceNotFoundException& snfe2) {
          //There was a first sequence, we just extend it:
          string ref1 = seq->getDescription();
          seq->setToSizeR(seq->size() + incomingBlock_->getNumberOfSites() + globalSpace);
          if (logstream_) {
            (*logstream_ << "BLOCK MERGER: extending " << ref1 << " with " << incomingBlock_->getNumberOfSites() << " gaps on the right.").endLine();
          }
        }
      } catch (SequenceNotFoundException& snfe1) {
        //There must be a second sequence then:
        seq.reset(new MafSequence(incomingBlock_->getSequenceForSpecies(allSp[i])));
        string ref2 = seq->getDescription();
        seq->setToSizeL(seq->size() + currentBlock_->getNumberOfSites() + globalSpace);
        if (logstream_) {
          (*logstream_ << "BLOCK MERGER: adding " << ref2 << " and extend it with " << currentBlock_->getNumberOfSites() << " gaps on the left.").endLine();
        }
      }
      mergedBlock->addSequence(*seq);
    }
    //Cleaning stuff:
    delete currentBlock_;
    delete incomingBlock_;
    currentBlock_ = mergedBlock;
    //We check if we can also merge the next block:
    incomingBlock_ = iterator_->nextBlock();
  }
  return currentBlock_;
}

