// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BlockMergerMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

std::unique_ptr<MafBlock> BlockMergerMafIterator::analyseCurrentBlock_()
{
  if (!incomingBlock_)
    return 0;
  currentBlock_  = std::move(incomingBlock_);
  incomingBlock_ = iterator_->nextBlock();
  while (incomingBlock_)
  {
    size_t globalSpace = 0;
    for (size_t i = 0; i < species_.size(); ++i)
    {
      try
      {
        const auto& seq1 = currentBlock_->sequenceForSpecies(species_[i]);
        const auto& seq2 = incomingBlock_->sequenceForSpecies(species_[i]);
        if (!seq1.hasCoordinates() || !seq2.hasCoordinates())
          throw Exception("BlockMergerMafIterator::nextBlock. Species '" + species_[i] + "' is missing coordinates in at least one block.");

        if (seq1.stop() > seq2.start())
          return std::move(currentBlock_);
        size_t space = seq2.start() - seq1.stop();
        if (space > maxDist_)
          return std::move(currentBlock_);
        if (i == 0)
          globalSpace = space;
        else
        {
          if (space != globalSpace)
            return std::move(currentBlock_);
        }
        if (seq1.getChromosome() != seq2.getChromosome()
            || VectorTools::contains(ignoreChrs_, seq1.getChromosome())
            || VectorTools::contains(ignoreChrs_, seq2.getChromosome())
            || seq1.getStrand() != seq2.getStrand()
            || seq1.getSrcSize() != seq2.getSrcSize())
        {
          // There is a syntheny break in this sequence, so we do not merge the blocks.
          return std::move(currentBlock_);
        }
      }
      catch (SequenceNotFoundException& snfe)
      {
        // At least one block does not contain the sequence.
        // We don't merge the blocks:
        return std::move(currentBlock_);
      }
    }
    // We merge the two blocks:
    if (logstream_)
    {
      (*logstream_ << "BLOCK MERGER: merging two consecutive blocks.").endLine();
    }
    vector<string> sp1 = currentBlock_->getSpeciesList();
    vector<string> sp2 = incomingBlock_->getSpeciesList();
    vector<string> allSp = VectorTools::unique(VectorTools::vectorUnion(sp1, sp2));
    // We need to create a new MafBlock:
    auto mergedBlock = make_unique<MafBlock>();
    // We average the score and pass values:
    unsigned int p1 = currentBlock_->getPass();
    unsigned int p2 = incomingBlock_->getPass();
    if (p1 == p2)
      mergedBlock->setPass(p1);
    double s1 = currentBlock_->getScore();
    double n1 = static_cast<double>(currentBlock_->getNumberOfSites());
    double s2 = incomingBlock_->getScore();
    double n2 = static_cast<double>(incomingBlock_->getNumberOfSites());
    mergedBlock->setScore((s1 * n1 + s2 * n2) / (n1 + n2));

    // Now fill the new block:
    for (size_t i = 0; i < allSp.size(); ++i)
    {
      unique_ptr<MafSequence> seq;
      try
      {
        seq = currentBlock_->removeSequenceForSpecies(allSp[i]);

        // Check is there is a second sequence:
        try
        {
          auto tmp = incomingBlock_->removeSequenceForSpecies(allSp[i]);
          string ref1 = seq->getDescription(), ref2 = tmp->getDescription();
          // Add spacer if needed:
          if (globalSpace > 0)
          {
            if (logstream_)
            {
              (*logstream_ << "BLOCK MERGER: a spacer of size " << globalSpace << " is inserted in sequence for species " << allSp[i] << ".").endLine();
            }
            seq->append(vector<int>(globalSpace, AlphabetTools::DNA_ALPHABET->getUnknownCharacterCode()));
          }
          if (seq->getChromosome() != tmp->getChromosome())
          {
            if (renameChimericChromosomes_)
            {
              if (seq->getChromosome().substr(0, 7) != "chimtig")
              {
                // Creates a new chimeric chromosome for this species:
                chimericChromosomeCounts_[seq->getSpecies()]++;
                seq->setChromosome("chimtig" + TextTools::toString(chimericChromosomeCounts_[seq->getSpecies()]));
              }
            }
            else
            {
              seq->setChromosome(seq->getChromosome() + "-" + tmp->getChromosome());
            }
            seq->removeCoordinates();
          }
          if (seq->getStrand() != tmp->getStrand())
          {
            seq->setStrand('?');
            seq->removeCoordinates();
          }
          if (seq->getName() != tmp->getName())
            tmp->setName(seq->getName()); // force name conversion to prevent exception in 'merge'.
          seq->merge(*tmp);
          if (logstream_)
          {
            (*logstream_ << "BLOCK MERGER: merging " << ref1 << " with " << ref2 << " into " << seq->getDescription()).endLine();
          }
        }
        catch (SequenceNotFoundException& snfe2)
        {
          // There was a first sequence, we just extend it:
          string ref1 = seq->getDescription();
          seq->setToSizeR(seq->size() + incomingBlock_->getNumberOfSites() + globalSpace);
          if (logstream_)
          {
            (*logstream_ << "BLOCK MERGER: extending " << ref1 << " with " << incomingBlock_->getNumberOfSites() << " gaps on the right.").endLine();
          }
        }
      }
      catch (SequenceNotFoundException& snfe1)
      {
        // There must be a second sequence then:
        seq = incomingBlock_->removeSequenceForSpecies(allSp[i]);
        string ref2 = seq->getDescription();
        seq->setToSizeL(seq->size() + currentBlock_->getNumberOfSites() + globalSpace);
        if (logstream_)
        {
          (*logstream_ << "BLOCK MERGER: adding " << ref2 << " and extend it with " << currentBlock_->getNumberOfSites() << " gaps on the left.").endLine();
        }
      }
      mergedBlock->addSequence(seq);
    }
    currentBlock_ = std::move(mergedBlock);
    // We check if we can also merge the next block:
    incomingBlock_ = iterator_->nextBlock();
  }
  return std::move(currentBlock_);
}
