//
// File: FeatureFilterMafIterator.cpp
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

#include "FeatureFilterMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* FeatureFilterMafIterator::analyseCurrentBlock_() throw (Exception)
{
  if (blockBuffer_.size() == 0) {
    //Unless there is no more block in the buffer, we need to parse more:
    do {
      MafBlock* block = iterator_->nextBlock();
      if (!block) return 0; //No more block.

      //Check if the block contains the reference species:
      if (!block->hasSequenceForSpecies(refSpecies_)) {
        if (logstream_) {
          (*logstream_ << "FEATURE FILTER: block " << block->getDescription() << " does not contain the reference species and was kept as is.").endLine(); 
        }
        return block;
      }

      //Get the feature ranges for this block:
      const MafSequence& refSeq = block->getSequenceForSpecies(refSpecies_);
      //first check if there is one (for now we assume that features refer to the chromosome or contig name, with implicit species):
      std::map<std::string, MultiRange<size_t> >::iterator mr = ranges_.find(refSeq.getChromosome());
      if (mr == ranges_.end()) {
        if (logstream_) {
          (*logstream_ << "FEATURE FILTER: block " << block->getDescription() << " does not contain any feature and was kept as is.").endLine(); 
        }
        return block;
      }
      //else
      MultiRange<size_t> mRange = mr->second;
      //mRange.restrictTo(Range<size_t>(refSeq.start(), refSeq.stop() + 1)); jdutheil on 17/04/13: do we really need the +1 here?
      mRange.restrictTo(refSeq.getRange(true));
      if (mRange.isEmpty()) {
        if (logstream_) {
          (*logstream_ << "FEATURE FILTER: block " << block->getDescription() << " does not contain any feature and was kept as is.").endLine(); 
        }
        return block;
      }
      std::vector<size_t> tmp = mRange.getBounds();

      //If the reference sequence is on the negative strand, then we have to correct the coordinates:
      std::deque<size_t> refBounds;
      if (refSeq.getStrand() == '-') {
        for (size_t i = 0; i < tmp.size(); ++i)
        {
          refBounds.push_front(refSeq.getSrcSize() - tmp[i]);
        }
      } else {
        refBounds = deque<size_t>(tmp.begin(), tmp.end());
      }

      //Now extract corresponding alignments. We use the range to split the original block.
      //Only thing to watch out is the coordinates, refering to the ref species...
      //A good idea is then to convert those with respect to the given block:

      int gap = refSeq.getAlphabet()->getGapCharacterCode();
      long int refPos = static_cast<long int>(refSeq.start()) - 1;
      //long int refPos = refSeq.getStrand() == '-' ? static_cast<long int>(refSeq.getSrcSize() - refSeq.start()) - 1 : static_cast<long int>(refSeq.start()) - 1;
      std::vector<size_t> pos;
      if (verbose_) {
        ApplicationTools::message->endLine();
        ApplicationTools::displayTask("Removing features", true);
      }
      for (size_t alnPos = 0; alnPos < refSeq.size() && refBounds.size() > 0; ++alnPos) {
        if (verbose_)
          ApplicationTools::displayGauge(alnPos, refBounds.back(), '>');
        if (refSeq[alnPos] != gap) {
          refPos++;
          //check if this position is a bound:
          if (refBounds.front() == static_cast<size_t>(refPos)) {
            pos.push_back(alnPos);
            refBounds.pop_front();
          }
        }
      }
      if (verbose_)
        ApplicationTools::displayTaskDone();

      //Check if the last bound matches the end of the alignment:
      if (refBounds.size() > 0 && refBounds.front() == refSeq.stop()) {
        pos.push_back(refSeq.size());
        refBounds.pop_front();
      }

      if (refBounds.size() > 0) {
        VectorTools::print(vector<size_t>(refBounds.begin(), refBounds.end()));
        throw Exception("FeatureFilterMafIterator::nextBlock(). An error occurred here, " + TextTools::toString(refBounds.size()) + " coordinates are left, in sequence " + refSeq.getDescription() + "... this is most likely a bug, please report!");
      }

      //Next step is simply to split the black according to the translated coordinates:
      if (pos.size() == 2 && pos.front() == 0 && pos.back() == block->getNumberOfSites()) {
        //Everything is removed:
        if (logstream_) {
          (*logstream_ << "FEATURE FILTER: block " << block->getDescription() << " was entirely removed. Tried to get the next one.").endLine();
        }
      } else {
        if (logstream_) {
          (*logstream_ << "FEATURE FILTER: block " << block->getDescription() << " with size "<< block->getNumberOfSites() << " will be split into " << (pos.size() / 2 + 1) << " blocks.").endLine();
        }
        if (verbose_) {
          ApplicationTools::displayTask("Spliting block", true);
        }
        for (size_t i = 0; i < pos.size(); i+=2) {
          if (verbose_)
            ApplicationTools::displayGauge(i, pos.size() - 2, '=');
          if (logstream_) {
            (*logstream_ << "FEATURE FILTER: removing region (" << pos[i] << ", " << pos[i+1] << ") from block " << block->getDescription() << ".").endLine();
          }
          if (pos[i] > 0) {
            MafBlock* newBlock = new MafBlock();
            newBlock->setScore(block->getScore());
            newBlock->setPass(block->getPass());
            for (size_t j = 0; j < block->getNumberOfSequences(); ++j) {
              MafSequence* subseq;
              if (i == 0) {
                subseq = block->getSequence(j).subSequence(0, pos[i]);
              } else {
                subseq = block->getSequence(j).subSequence(pos[i - 1], pos[i] - pos[i - 1]);
              }
              newBlock->addSequence(*subseq);
              delete subseq;
            }
            blockBuffer_.push_back(newBlock);
          }
        
          if (keepTrashedBlocks_) {
            MafBlock* outBlock = new MafBlock();
            outBlock->setScore(block->getScore());
            outBlock->setPass(block->getPass());
            for (size_t j = 0; j < block->getNumberOfSequences(); ++j) {
              MafSequence* outseq = block->getSequence(j).subSequence(pos[i], pos[i + 1] - pos[i]);
              outBlock->addSequence(*outseq);
              delete outseq;
            } 
            trashBuffer_.push_back(outBlock);
          }
        }
        //Add last block:
        if (pos.back() < block->getNumberOfSites()) {
          MafBlock* newBlock = new MafBlock();
          newBlock->setScore(block->getScore());
          newBlock->setPass(block->getPass());
          for (size_t j = 0; j < block->getNumberOfSequences(); ++j) {
            MafSequence* subseq;
            subseq = block->getSequence(j).subSequence(pos[pos.size() - 1], block->getNumberOfSites() - pos[pos.size() - 1]);
            newBlock->addSequence(*subseq);
            delete subseq;
          }
          blockBuffer_.push_back(newBlock);
        }
        if (verbose_)
          ApplicationTools::displayTaskDone();

        delete block;
      }
    } while (blockBuffer_.size() == 0);
  }

  MafBlock* nxtBlock = blockBuffer_.front();
  blockBuffer_.pop_front();
  return nxtBlock;
}

