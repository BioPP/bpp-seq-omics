// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ConcatenateMafIterator.h"

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

// From bpp-core:
#include <Bpp/App/ApplicationTools.h>

using namespace std;

unique_ptr<MafBlock> ConcatenateMafIterator::analyseCurrentBlock_()
{
  if (!incomingBlock_)
    return 0;
  currentBlock_  = std::move(incomingBlock_);
  incomingBlock_ = iterator_->nextBlock();
  size_t count = 1;
  if (verbose_)
    ApplicationTools::displayMessage("Concatenating new block...");
  while (incomingBlock_ &&
      (refSpecies_ == "" ||
      (incomingBlock_->hasSequenceForSpecies(refSpecies_) &&
      currentBlock_->hasSequenceForSpecies(refSpecies_) &&
      incomingBlock_->sequenceForSpecies(refSpecies_).getChromosome() ==
      currentBlock_->sequenceForSpecies(refSpecies_).getChromosome()
      )
      )
      )
  {
    if (currentBlock_->getNumberOfSites() >= minimumSize_)
    {
      return std::move(currentBlock_);
    }
    if (verbose_)
    {
      ApplicationTools::displayUnlimitedGauge(count++, "Concatenating...");
    }

    // We merge the two blocks:
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
        seq.reset(new MafSequence(currentBlock_->sequenceForSpecies(allSp[i])));

        // Check is there is a second sequence:
        try
        {
          auto tmp = make_unique<MafSequence>(incomingBlock_->sequenceForSpecies(allSp[i]));
          string ref1 = seq->getDescription(), ref2 = tmp->getDescription();
          if (seq->getChromosome() != tmp->getChromosome())
          {
            seq->setChromosome("fus");
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
            (*logstream_ << "BLOCK CONCATENATE: merging " << ref1 << " with " << ref2 << " into " << seq->getDescription()).endLine();
          }
        }
        catch (SequenceNotFoundException& snfe2)
        {
          // There was a first sequence, we just extend it:
          string ref1 = seq->getDescription();
          seq->setToSizeR(seq->size() + incomingBlock_->getNumberOfSites());
          if (logstream_)
          {
            (*logstream_ << "BLOCK CONCATENATE: extending " << ref1 << " with " << incomingBlock_->getNumberOfSites() << " gaps on the right.").endLine();
          }
        }
      }
      catch (SequenceNotFoundException& snfe1)
      {
        // There must be a second sequence then:
        seq.reset(new MafSequence(incomingBlock_->sequenceForSpecies(allSp[i])));
        string ref2 = seq->getDescription();
        seq->setToSizeL(seq->size() + currentBlock_->getNumberOfSites());
        if (logstream_)
        {
          (*logstream_ << "BLOCK CONCATENATE: adding " << ref2 << " and extend it with " << currentBlock_->getNumberOfSites() << " gaps on the left.").endLine();
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
