//
// File: FeatureExtractorMafIterator.cpp
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

#include "FeatureExtractorMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* FeatureExtractorMafIterator::analyseCurrentBlock_() throw (Exception)
{
  while (blockBuffer_.size() == 0) {
    //Unless there is no more block in the buffer, we need to parse more:
    unique_ptr<MafBlock> block;
    START:
    block.reset(iterator_->nextBlock());
    if (!block.get()) return 0; //No more block.

    //Check if the block contains the reference species:
    if (!block->hasSequenceForSpecies(refSpecies_))
      goto START;

    //Get the feature ranges for this block:
    const MafSequence& refSeq = block->getSequenceForSpecies(refSpecies_);
    //first check if there is one (for now we assume that features refer to the chromosome or contig name, with implicit species):
    std::map<std::string, RangeSet<size_t> >::iterator mr = ranges_.find(refSeq.getChromosome());
    if (mr == ranges_.end())
      goto START;
        
    RangeSet<size_t> ranges = mr->second;
    if (completeOnly_)
      ranges.filterWithin(refSeq.getRange(true));
    else  
      ranges.restrictTo(refSeq.getRange(true));
    if (ranges.isEmpty())
      goto START;

    //If the reference sequence is on the negative strand, then we have to correct the coordinates:
    (*logstream_ << "Strand: " << refSeq.getStrand()).endLine();
    (*logstream_ << refSeq.start()).endLine();
    (*logstream_ << "ok :)").endLine();
    if (refSeq.getStrand() == '-') {
      RangeSet<size_t> cRanges;
      for (set<Range<size_t>*>::iterator it = ranges.getSet().begin();
          it != ranges.getSet().end();
          ++it)
      {
        cRanges.addRange(SeqRange(refSeq.getSrcSize() - (**it).end(), refSeq.getSrcSize() - (**it).begin(), dynamic_cast<SeqRange*>(*it)->getStrand()));
      }
      ranges = cRanges;
    }

    //We will need to convert to alignment positions, using a sequence walker:
    SequenceWalker walker(refSeq);

    //Now creates all blocks for all ranges:
    if (verbose_) {
      ApplicationTools::message->endLine();
      ApplicationTools::displayTask("Extracting annotations", true);
    }
    if (logstream_) {
      (*logstream_ << "FEATURE EXTRACTOR: extracting " << ranges.getSet().size() << " features from block " << block->getDescription() << ".").endLine();
    }

    size_t i = 0;
    for (set<Range<size_t>*>::iterator it = ranges.getSet().begin();
        it !=  ranges.getSet().end();
        ++it)
    {
      if (verbose_) {
        ApplicationTools::displayGauge(i++, ranges.getSet().size() - 1, '=');
      } 
      //This does not go after i=0, problem with ranges?????
      MafBlock* newBlock = new MafBlock();
      newBlock->setScore(block->getScore());
      newBlock->setPass(block->getPass());
      size_t a = walker.getAlignmentPosition((**it).begin() - refSeq.start());
      size_t b = walker.getAlignmentPosition((**it).end() - refSeq.start() - 1);
      for (size_t j = 0; j < block->getNumberOfSequences(); ++j) {
        unique_ptr<MafSequence> subseq;
        subseq.reset(block->getSequence(j).subSequence(a, b - a + 1));
        if (!ignoreStrand_) {
          if ((dynamic_cast<SeqRange*>(*it)->isNegativeStrand() && refSeq.getStrand() == '+') ||
             (!dynamic_cast<SeqRange*>(*it)->isNegativeStrand() && refSeq.getStrand() == '-'))
          {
            SequenceTools::invertComplement(*subseq);
          }
        }
        newBlock->addSequence(*subseq);
      }
      blockBuffer_.push_back(newBlock);
    }
        
    if (verbose_)
      ApplicationTools::displayTaskDone();

  }

  MafBlock* nxtBlock = blockBuffer_.front();
  blockBuffer_.pop_front();
  return nxtBlock;
}

