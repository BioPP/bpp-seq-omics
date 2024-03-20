// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "EntropyFilterMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> EntropyFilterMafIterator::analyseCurrentBlock_()
{
  if (blockBuffer_.size() == 0)
  {
    // Else there is no more block in the buffer, we need to parse more:
    do
    {
      auto block = iterator_->nextBlock();
      if (!block)
        return 0; // No more block.

      // Parse block.
      int gap = AlphabetTools::DNA_ALPHABET->getGapCharacterCode();
      // unused for now
      // int unk = AlphabetTools::DNA_ALPHABET.getUnknownCharacterCode();
      size_t nr;
      size_t nc = static_cast<size_t>(block->getNumberOfSites());
      if (nc < windowSize_)
        throw Exception("EntropyFilterMafIterator::analyseCurrentBlock_. Block is smaller than window size: " + TextTools::toString(nc));

      vector< vector<int>> aln;
      if (missingAsGap_ && !ignoreGaps_)
      {
        nr = species_.size();
        aln.resize(nr);
        for (size_t i = 0; i < nr; ++i)
        {
          if (block->hasSequenceForSpecies(species_[i]))
            aln[i] = block->sequenceForSpecies(species_[i]).getContent();
          else
          {
            aln[i].resize(nc);
            fill(aln[i].begin(), aln[i].end(), gap);
          }
        }
      }
      else
      {
        vector<string> speciesSet = VectorTools::vectorIntersection(species_, block->getSpeciesList());
        nr = speciesSet.size();
        aln.resize(nr);
        for (size_t i = 0; i < nr; ++i)
        {
          aln[i] = block->sequenceForSpecies(species_[i]).getContent();
        }
      }
      // First we create a mask:
      vector<size_t> pos;
      // Reset window:
      window_.clear();
      // Init window:
      size_t i;
      for (i = 0; i < windowSize_; ++i)
      {
        vector<int> col;
        for (size_t j = 0; j < nr; ++j)
        {
          int x = aln[j][i];
          if (x != gap || !ignoreGaps_)
            col.push_back(x);
        }
        double entropy = VectorTools::shannonDiscrete<int, double>(col) / log(5.);
        window_.push_back(entropy > maxEnt_ ? 1 : 0);
      }
      // Slide window:
      if (verbose_)
      {
        ApplicationTools::message->endLine();
        ApplicationTools::displayTask("Sliding window for entropy filter", true);
      }
      while (i + step_ < nc)
      {
        if (verbose_)
          ApplicationTools::displayGauge(i - windowSize_, nc - windowSize_ - 1, '>');
        // Evaluate current window:
        unsigned int count = std::accumulate(window_.begin(), window_.end(), 0u);
        if (count > maxPos_)
        {
          if (pos.size() == 0)
          {
            pos.push_back(i - windowSize_);
            pos.push_back(i);
          }
          else
          {
            if (i - windowSize_ < pos[pos.size() - 1])
            {
              pos[pos.size() - 1] = i; // Windows are overlapping and we extend previous region
            }
            else // This is a new region
            {
              pos.push_back(i - windowSize_);
              pos.push_back(i);
            }
          }
        }

        // Move forward:
        for (size_t k = 0; k < step_; ++k)
        {
          vector<int> col;
          for (size_t j = 0; j < nr; ++j)
          {
            int x = aln[j][i];
            if (x != gap || !ignoreGaps_)
              col.push_back(x);
          }
          double entropy = VectorTools::shannonDiscrete<int, double>(col) / log(5.);
          window_.push_back(entropy > maxEnt_ ? 1 : 0);
          window_.pop_front();
          ++i;
        }
      }

      // Evaluate last window:
      unsigned int count = std::accumulate(window_.begin(), window_.end(), 0u);
      if (count > maxPos_)
      {
        if (pos.size() == 0)
        {
          pos.push_back(i - windowSize_);
          pos.push_back(i);
        }
        else
        {
          if (i - windowSize_ <= pos[pos.size() - 1])
          {
            pos[pos.size() - 1] = i; // Windows are overlapping and we extend previous region
          }
          else // This is a new region
          {
            pos.push_back(i - windowSize_);
            pos.push_back(i);
          }
        }
      }
      if (verbose_)
        ApplicationTools::displayTaskDone();

      // Now we remove regions with two many gaps, using a sliding window:
      if (pos.size() == 0)
      {
        blockBuffer_.push_back(std::move(block));
        if (logstream_)
        {
          (*logstream_ << "ENTROPY CLEANER: block " << block->getDescription() << " is clean and kept as is.").endLine();
        }
      }
      else if (pos.size() == 2 && pos.front() == 0 && pos.back() == block->getNumberOfSites())
      {
        // Everything is removed:
        if (logstream_)
        {
          (*logstream_ << "ENTROPY CLEANER: block " << block->getDescription() << " was entirely removed. Tried to get the next one.").endLine();
        }
      }
      else
      {
        if (logstream_)
        {
          (*logstream_ << "ALN CLEANER: block " << block->getDescription() << " with size " << block->getNumberOfSites() << " will be split into " << (pos.size() / 2 + 1) << " blocks.").endLine();
        }
        if (verbose_)
        {
          ApplicationTools::message->endLine();
          ApplicationTools::displayTask("Spliting block", true);
        }
        for (i = 0; i < pos.size(); i += 2)
        {
          if (verbose_)
            ApplicationTools::displayGauge(i, pos.size() - 2, '=');
          if (logstream_)
          {
            (*logstream_ << "ENTROPY CLEANER: removing region (" << pos[i] << ", " << pos[i + 1] << ") from block " << block->getDescription() << ".").endLine();
          }
          if (pos[i] > 0)
          {
            auto newBlock = make_unique<MafBlock>();
            newBlock->setScore(block->getScore());
            newBlock->setPass(block->getPass());
            for (size_t j = 0; j < block->getNumberOfSequences(); ++j)
            {
              unique_ptr<MafSequence> subseq;
              if (i == 0)
              {
                subseq = block->sequence(j).subSequence(0, pos[i]);
              }
              else
              {
                subseq = block->sequence(j).subSequence(pos[i - 1], pos[i] - pos[i - 1]);
              }
              newBlock->addSequence(subseq);
            }
            blockBuffer_.push_back(std::move(newBlock));
          }

          if (keepTrashedBlocks_)
          {
            auto outBlock = make_unique<MafBlock>();
            outBlock->setScore(block->getScore());
            outBlock->setPass(block->getPass());
            for (size_t j = 0; j < block->getNumberOfSequences(); ++j)
            {
              auto outseq = block->sequence(j).subSequence(pos[i], pos[i + 1] - pos[i]);
              outBlock->addSequence(outseq);
            }
            trashBuffer_.push_back(std::move(outBlock));
          }
        }
        // Add last block:
        if (pos[pos.size() - 1] < block->getNumberOfSites())
        {
          auto newBlock = make_unique<MafBlock>();
          newBlock->setScore(block->getScore());
          newBlock->setPass(block->getPass());
          for (size_t j = 0; j < block->getNumberOfSequences(); ++j)
          {
            auto subseq = block->sequence(j).subSequence(pos[pos.size() - 1], block->getNumberOfSites() - pos[pos.size() - 1]);
            newBlock->addSequence(subseq);
          }
          blockBuffer_.push_back(std::move(newBlock));
        }
        if (verbose_)
          ApplicationTools::displayTaskDone();
      }
    }
    while (blockBuffer_.size() == 0);
  }

  auto block = std::move(blockBuffer_.front());
  blockBuffer_.pop_front();
  return block;
}
