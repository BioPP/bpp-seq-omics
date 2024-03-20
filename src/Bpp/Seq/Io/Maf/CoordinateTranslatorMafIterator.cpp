// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CoordinateTranslatorMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> CoordinateTranslatorMafIterator::analyseCurrentBlock_()
{
  auto block = iterator_->nextBlock();
  if (!block)
    return nullptr; // No more block.

  // Check if the block contains the reference and target species:
  if (!block->hasSequenceForSpecies(referenceSpecies_))
    return block;
  if (!block->hasSequenceForSpecies(targetSpecies_))
    return block;

  // Get the feature ranges for this block:
  const auto& refSeq = block->sequenceForSpecies(referenceSpecies_);
  const auto& targetSeq = block->sequenceForSpecies(targetSpecies_);

  // first check if there is one (for now we assume that features refer to the chromosome or contig name, with implicit species):
  std::map<std::string, SequenceFeatureSet*>::iterator mr = inputFeaturesPerChr_.find(refSeq.getChromosome());
  if (mr == inputFeaturesPerChr_.end())
    return block;

  // second get only features within this block:
  unique_ptr<SequenceFeatureSet> selectedFeatures(mr->second->getSubsetForRange(SeqRange(refSeq.getRange(true)), true));

  // test if there are some features to translate here:
  if (selectedFeatures->isEmpty())
    return block;

  // Get coordinate range sets:
  RangeSet<size_t> ranges;
  selectedFeatures->fillRangeCollection(ranges);

  // If the reference sequence is on the negative strand, then we have to correct the coordinates:
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
  SequenceWalker referenceWalker(refSeq);
  SequenceWalker targetWalker(targetSeq);

  // Now creates all blocks for all ranges:
  if (verbose_)
  {
    ApplicationTools::message->endLine();
    ApplicationTools::displayTask("Extracting annotations", true);
  }
  if (logstream_)
  {
    (*logstream_ << "COORDINATE CONVERTOR: lifting over " << ranges.getSet().size() << " features from block " << block->getDescription() << ".").endLine();
  }

  auto alphabet = refSeq.getAlphabet();
  size_t i = 0;
  for (const auto& it : ranges.getSet())
  {
    if (verbose_)
    {
      ApplicationTools::displayGauge(i++, ranges.getSet().size() - 1, '=');
    }
    size_t a = referenceWalker.getAlignmentPosition(it->begin() - refSeq.start());
    size_t b = referenceWalker.getAlignmentPosition(it->end() - refSeq.start() - 1);
    string targetPos1 = "NA", targetPos2 = "NA";
    if (!alphabet->isGap(targetSeq[a]) || outputClosestCoordinate_)
    {
      size_t a2 = targetWalker.getSequencePosition(a) + targetSeq.start();
      if (targetSeq.getStrand() == '-')
      {
        a2 = targetSeq.getSrcSize() - a2;
      }
      targetPos1 = TextTools::toString(a2);
    }
    if (!alphabet->isGap(targetSeq[b]) || outputClosestCoordinate_)
    {
      size_t b2 = targetWalker.getSequencePosition(b) + targetSeq.start() + 1;
      if (targetSeq.getStrand() == '-')
      {
        b2 = targetSeq.getSrcSize() - b2;
      }
      targetPos2 = TextTools::toString(b2);
    }
    output_ << refSeq.getChromosome() << "\t" << refSeq.getStrand() << "\t" << it->begin() << "\t" << it->end() << "\t";
    output_ << targetSeq.getChromosome() << "\t" << targetSeq.getStrand() << "\t" << targetPos1 << "\t" << targetPos2 << endl;
  }

  if (verbose_)
    ApplicationTools::displayTaskDone();

  // Block is simply forwarded:
  return block;
}
