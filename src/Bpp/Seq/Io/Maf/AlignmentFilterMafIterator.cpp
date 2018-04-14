//
// File: AlignmentFilterMafIterator.cpp
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

#include "AlignmentFilterMafIterator.h"

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* AlignmentFilterMafIterator::analyseCurrentBlock_()
{
  if (blockBuffer_.size() == 0) {
    //Else there is no more block in the buffer, we need to parse more:
    do {
      MafBlock* block = iterator_->nextBlock();
      if (!block) return 0; //No more block.
    
      //Parse block.
      int gap = AlphabetTools::DNA_ALPHABET.getGapCharacterCode();
      int unk = AlphabetTools::DNA_ALPHABET.getUnknownCharacterCode();
      size_t nr;
      size_t nc = static_cast<size_t>(block->getNumberOfSites());
      if (nc < windowSize_)
        throw Exception("AlignmentFilterMafIterator::analyseCurrentBlock_. Block is smaller than window size: " + TextTools::toString(nc));

      vector< vector<int> > aln;
      if (missingAsGap_) {
        nr = species_.size();
        aln.resize(nr);
        for (size_t i = 0; i < nr; ++i) {
          if (block->hasSequenceForSpecies(species_[i]))
            aln[i] = block->getSequenceForSpecies(species_[i]).getContent();
          else {
            aln[i].resize(nc); 
            fill(aln[i].begin(), aln[i].end(), gap);
          } 
        }
      } else {
        vector<string> speciesSet = VectorTools::vectorIntersection(species_, block->getSpeciesList());
        for (size_t i = 0; i < species_.size(); ++i) {
          if (block->hasSequenceForSpecies(species_[i])) {
            aln.push_back(block->getSequenceForSpecies(species_[i]).getContent());
          } else {
            if (!relative_) {
              throw Exception("AlignmentFilterMafIterator::analyseCurrentBlock_. Block does not include selected species '" + species_[i] + "' and threshold are absolutes, leading to an undefined behavior. Consider selecting blocks first, or use relative thresholds.");
            }
          }
        }
        nr = aln.size();
      }

      //First we create a mask:
      vector<size_t> pos;
      vector<int> col(nr);
      //Reset window:
      window_.clear();
      //Init window:
      size_t i;
      for (i = 0; i < windowSize_; ++i) {
        for (size_t j = 0; j < nr; ++j) {
          col[j] = aln[j][i];
        }
        window_.push_back(col);
      }
      //Slide window:
      if (verbose_) {
        ApplicationTools::message->endLine();
        ApplicationTools::displayTask("Sliding window for alignment filter", true);
      }
      while (i + step_ < nc) {
        if (verbose_)
          ApplicationTools::displayGauge(i - windowSize_, nc - windowSize_ - 1, '>');
        //Evaluate current window:
        unsigned int sumGap = 0;
        double sumEnt = 0;
        for (size_t u = 0; u < window_.size(); ++u) {
          for (size_t v = 0; v < window_[u].size(); ++v) {
            if (window_[u][v] == gap || window_[u][v] == unk) {
              sumGap++;
            }
            col[v] = window_[u][v];
            if (col[v] == unk) col[v] = gap;
          }
          sumEnt += VectorTools::shannonDiscrete<int, double>(col) / log(5.);
        }
        bool test = (sumEnt / static_cast<double>(windowSize_)) > maxEnt_;
        if (relative_) {
          double propGap = static_cast<double>(sumGap) / static_cast<double>(nr * windowSize_);
          test = test && (propGap > maxPropGap_);
        } else {
          test = test && (sumGap > maxGap_);
        }
        if (test) {
          if (pos.size() == 0) {
            pos.push_back(i - windowSize_);
            pos.push_back(i);
          } else {
            if (i - windowSize_ <= pos[pos.size() - 1]) {
              pos[pos.size() - 1] = i; //Windows are overlapping and we extend previous region
            } else { //This is a new region
              pos.push_back(i - windowSize_);
              pos.push_back(i);
            }
          }
        }
      
        //Move forward:
        for (size_t k = 0; k < step_; ++k) {
          for (size_t j = 0; j < nr; ++j) {
            col[j] = aln[j][i];
          }
          window_.push_back(col);
          window_.pop_front();
          ++i;
        }
      }
      
      //Evaluate last window:
      unsigned int sumGap = 0;
      double sumEnt = 0;
      for (size_t u = 0; u < window_.size(); ++u) {
        for (size_t v = 0; v < window_[u].size(); ++v) {
          if (window_[u][v] == gap || window_[u][v] == unk) {
            sumGap++;
          }
          col[v] = window_[u][v];
          if (col[v] == unk) col[v] = gap;
        }
        sumEnt += VectorTools::shannonDiscrete<int, double>(col) / log(5.);
      }
      bool test = (sumEnt / static_cast<double>(windowSize_)) > maxEnt_;
      if (relative_) {
        double propGap = static_cast<double>(sumGap) / static_cast<double>(nr * windowSize_);
        test = test && (propGap > maxPropGap_);
      } else {
        test = test && (sumGap > maxGap_);
      }
      if (test) {
        if (pos.size() == 0) {
          pos.push_back(i - windowSize_);
          pos.push_back(i);
        } else {
          if (i - windowSize_ <= pos[pos.size() - 1]) {
            pos[pos.size() - 1] = i; //Windows are overlapping and we extend previous region
          } else { //This is a new region
            pos.push_back(i - windowSize_);
            pos.push_back(i);
          }
        }
      } 
      if (verbose_)
        ApplicationTools::displayTaskDone();
    
      //Now we remove regions with two many gaps, using a sliding window:
      if (pos.size() == 0) {
        blockBuffer_.push_back(block);
        if (logstream_) {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " is clean and kept as is.").endLine();
        }
      } else if (pos.size() == 2 && pos.front() == 0 && pos.back() == block->getNumberOfSites()) {
        //Everything is removed:
        if (logstream_) {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " was entirely removed. Tried to get the next one.").endLine();
        }
      } else {
        if (logstream_) {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " with size "<< block->getNumberOfSites() << " will be split into " << (pos.size() / 2 + 1) << " blocks.").endLine();
        }
        if (verbose_) {
          ApplicationTools::message->endLine();
          ApplicationTools::displayTask("Spliting block", true);
        }
        for (i = 0; i < pos.size(); i+=2) {
          if (verbose_)
            ApplicationTools::displayGauge(i, pos.size() - 2, '=');
          if (logstream_) {
            (*logstream_ << "ALN CLEANER: removing region (" << pos[i] << ", " << pos[i+1] << ") from block " << block->getDescription() << ".").endLine();
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
        if (pos[pos.size() - 1] < block->getNumberOfSites()) {
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

  MafBlock* block = blockBuffer_.front();
  blockBuffer_.pop_front();
  return block;
}

MafBlock* AlignmentFilter2MafIterator::analyseCurrentBlock_()
{
  if (blockBuffer_.size() == 0) {
    //Else there is no more block in the buffer, we need to parse more:
    do {
      MafBlock* block = iterator_->nextBlock();
      if (!block) return 0; //No more block.
    
      //Parse block.
      int gap = AlphabetTools::DNA_ALPHABET.getGapCharacterCode();
      int unk = AlphabetTools::DNA_ALPHABET.getUnknownCharacterCode();
      size_t nr;
      size_t nc = static_cast<size_t>(block->getNumberOfSites());
      if (nc < windowSize_)
        throw Exception("AlignmentFilter2MafIterator::analyseCurrentBlock_. Block is smaller than window size: " + TextTools::toString(nc));

      vector< vector<int> > aln;
      if (missingAsGap_) {
        nr = species_.size();
        aln.resize(nr);
        for (size_t i = 0; i < nr; ++i) {
          if (block->hasSequenceForSpecies(species_[i]))
            aln[i] = block->getSequenceForSpecies(species_[i]).getContent();
          else {
            aln[i].resize(nc); 
            fill(aln[i].begin(), aln[i].end(), gap);
          } 
        }
      } else {
        vector<string> speciesSet = VectorTools::vectorIntersection(species_, block->getSpeciesList());
        for (size_t i = 0; i < species_.size(); ++i) {
          if (block->hasSequenceForSpecies(species_[i])) {
            aln.push_back(block->getSequenceForSpecies(species_[i]).getContent());
          } else {
            if (!relative_) {
              throw Exception("AlignmentFilter2MafIterator::analyseCurrentBlock_. Block does not include selected species '" + species_[i] + "' and threshold are absolutes, leading to an undefined behavior. Consider selecting blocks first, or use relative thresholds.");
            }
          }
        }
        nr = aln.size();
      }

      //First we create a mask:
      vector<size_t> pos;
      vector<bool> col(nr);
      //Reset window:
      window_.clear();
      //Init window:
      size_t i;
      for (i = 0; i < windowSize_; ++i) {
        for (size_t j = 0; j < nr; ++j) {
          col[j] = (aln[j][i] == gap || aln[j][i] == unk);
        }
        window_.push_back(col);
      }
      //Slide window:
      if (verbose_) {
        ApplicationTools::message->endLine();
        ApplicationTools::displayTask("Sliding window for alignment filter", true);
      }
      while (i + step_ < nc) {
        if (verbose_)
          ApplicationTools::displayGauge(i - windowSize_, nc - windowSize_ - 1, '>');
        //Evaluate current window:
        unsigned int count = 0;
        bool posIsGap = false;
        for (size_t u = 0; u < window_.size(); ++u) {
          unsigned int partialCount = 0;
          if (!posIsGap || (u > 0 && window_[u] != window_[u - 1])) {
            for (size_t v = 0; v < window_[u].size(); ++v)
              if (window_[u][v]) partialCount++;
            bool test;
            if (relative_) {
              test = (static_cast<double>(partialCount) / static_cast<double>(nr) > maxPropGap_);
            } else {
              test = (partialCount > maxGap_);
            }
            if (test) {
              count++;
              posIsGap = true;
            } else {
              posIsGap = false;
            }
          }
        }
        if (count > maxPos_) {
          if (pos.size() == 0) {
            pos.push_back(i - windowSize_);
            pos.push_back(i);
          } else {
            if (i - windowSize_ <= pos[pos.size() - 1]) {
              pos[pos.size() - 1] = i; //Windows are overlapping and we extend previous region
            } else { //This is a new region
              pos.push_back(i - windowSize_);
              pos.push_back(i);
            }
          }
        }
      
        //Move forward:
        for (size_t k = 0; k < step_; ++k) {
          for (size_t j = 0; j < nr; ++j) {
            col[j] = (aln[j][i] == gap || aln[j][i] == unk);
          }
          window_.push_back(col);
          window_.pop_front();
          ++i;
        }
      }

      //Evaluate last window:
      unsigned int count = 0;
      bool posIsGap = false;
      for (size_t u = 0; u < window_.size(); ++u) {
        unsigned int partialCount = 0;
        if (!posIsGap || (u > 0 && window_[u] != window_[u - 1])) {
          for (size_t v = 0; v < window_[u].size(); ++v)
            if (window_[u][v]) partialCount++;
          bool test;
          if (relative_) {
            test = (static_cast<double>(partialCount) / static_cast<double>(nr) > maxPropGap_);
          } else {
            test = (partialCount > maxGap_);
          }
          if (test) {
            count++;
            posIsGap = true;
          } else {
            posIsGap = false;
          }
        }
      }
      if (count > maxPos_) {
        if (pos.size() == 0) {
          pos.push_back(i - windowSize_);
          pos.push_back(i);
        } else {
          if (i - windowSize_ <= pos[pos.size() - 1]) {
            pos[pos.size() - 1] = i; //Windows are overlapping and we extend previous region
          } else { //This is a new region
            pos.push_back(i - windowSize_);
            pos.push_back(i);
          }
        }
      } 
      if (verbose_)
        ApplicationTools::displayTaskDone();
    
      //Now we remove regions with two many gaps, using a sliding window:
      if (pos.size() == 0) {
        blockBuffer_.push_back(block);
        if (logstream_) {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " is clean and kept as is.").endLine();
        }
      } else if (pos.size() == 2 && pos.front() == 0 && pos.back() == block->getNumberOfSites()) {
        //Everything is removed:
        if (logstream_) {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " was entirely removed. Tried to get the next one.").endLine();
        }
      } else {
        if (logstream_) {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " with size "<< block->getNumberOfSites() << " will be split into " << (pos.size() / 2 + 1) << " blocks.").endLine();
        }
        if (verbose_) {
          ApplicationTools::message->endLine();
          ApplicationTools::displayTask("Spliting block", true);
        }
        for (i = 0; i < pos.size(); i+=2) {
          if (verbose_)
            ApplicationTools::displayGauge(i, pos.size() - 2, '=');
          if (logstream_) {
            (*logstream_ << "ALN CLEANER: removing region (" << pos[i] << ", " << pos[i+1] << ") from block " << block->getDescription() << ".").endLine();
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
        if (pos[pos.size() - 1] < block->getNumberOfSites()) {
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

  MafBlock* block = blockBuffer_.front();
  blockBuffer_.pop_front();
  return block;
}

