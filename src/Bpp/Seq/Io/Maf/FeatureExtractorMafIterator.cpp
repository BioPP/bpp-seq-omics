// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "FeatureExtractorMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> FeatureExtractorMafIterator::analyseCurrentBlock_()
{
  while (blockBuffer_.size() == 0)
  {
    // Unless there is no more block in the buffer, we need to parse more:
    unique_ptr<MafBlock> block;
START:
    block = iterator_->nextBlock();
    if (!block.get())
      return 0; // No more block.

    // Check if the block contains the reference species:
    if (!block->hasSequenceForSpecies(refSpecies_)) {
      goto START;
    }

    // Get the feature ranges for this block:
    const auto& refSeq = block->sequenceForSpecies(refSpecies_);
    // first check if there is one (for now we assume that features refer to the chromosome or contig name, with implicit species):

    auto mr = ranges_.find(refSeq.getChromosome());
    if (mr == ranges_.end()) {
      goto START;
    }

    RangeSet<size_t> ranges = mr->second;
    if (completeOnly_)
      ranges.filterWithin(refSeq.getRange(true));
    else
      ranges.restrictTo(refSeq.getRange(true));
    if (ranges.isEmpty()) {
      goto START;
    }

    // If the reference sequence is on the negative strand, then we have to correct the coordinates:
    (*logstream_ << "Strand: " << refSeq.getStrand()).endLine();
    if (refSeq.getStrand() == '-')
    {
      RangeSet<size_t> cRanges;
      for (const auto& it : ranges.getSet())
      {
        cRanges.addRange(SeqRange(refSeq.getSrcSize() - it->end(), refSeq.getSrcSize() - it->begin(), dynamic_cast<const SeqRange*>(it)->getStrand()));
      }
      ranges = cRanges;
    }

    // We will need to convert to alignment positions, using a sequence walker:
    SequenceWalker walker(refSeq);

    // Now creates all blocks for all ranges:
    if (verbose_)
    {
      ApplicationTools::message->endLine();
      ApplicationTools::displayTask("Extracting annotations", true);
    }
    if (logstream_)
    {
      (*logstream_ << "FEATURE EXTRACTOR: extracting " << ranges.getSet().size() << " features from block " << block->getDescription() << ".").endLine();
    }

    size_t i = 0;
    for (const auto& it : ranges.getSet())
    {
      if (verbose_)
      {
        ApplicationTools::displayGauge(i++, ranges.getSet().size() - 1, '=');
      }
      // This does not go after i=0, problem with ranges?????
      auto newBlock = make_unique<MafBlock>();
      newBlock->setScore(block->getScore());
      newBlock->setPass(block->getPass());
      size_t a = walker.getAlignmentPosition(it->begin() - refSeq.start());
      size_t b = walker.getAlignmentPosition(it->end() - refSeq.start() - 1);
      for (size_t j = 0; j < block->getNumberOfSequences(); ++j)
      {
        auto subseq = block->sequence(j).subSequence(a, b - a + 1);
        if (!ignoreStrand_)
        {
          if ((dynamic_cast<const SeqRange*>(it)->isNegativeStrand() && refSeq.getStrand() == '+') ||
              (!dynamic_cast<const SeqRange*>(it)->isNegativeStrand() && refSeq.getStrand() == '-'))
          {
            SequenceTools::invertComplement(*subseq);
          }
        }
        (*logstream_ << subseq->getName()).endLine();
        newBlock->addSequence(subseq);
      }
      blockBuffer_.push_back(move(newBlock));
    }

    if (verbose_)
      ApplicationTools::displayTaskDone();
  }

  auto nxtBlock = move(blockBuffer_.front());
  blockBuffer_.pop_front();
  return nxtBlock;
}
