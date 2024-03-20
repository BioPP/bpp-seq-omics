// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OutputAlignmentMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> OutputAlignmentMafIterator::analyseCurrentBlock_()
{
  auto block = iterator_->nextBlock();
  if (block)
  {
    if (output_)
    {
      writeBlock(*output_, *block);
    }
    else
    {
      string chr   = "ChrNA";
      string start = "StartNA";
      string stop  = "StopNA";
      if (block->hasSequenceForSpecies(refSpecies_))
      {
        const auto& refseq = block->sequenceForSpecies(refSpecies_);
        chr   = refseq.getChromosome();
        start = TextTools::toString(refseq.start());
        stop  = TextTools::toString(refseq.stop());
      }
      string file = file_;
      TextTools::replaceAll(file, "%i", TextTools::toString(++currentBlockIndex_));
      TextTools::replaceAll(file, "%c", chr);
      TextTools::replaceAll(file, "%b", start);
      TextTools::replaceAll(file, "%e", stop);
      std::ofstream output(file.c_str(), ios::out);

      writeBlock(output, *block);
    }
  }
  return block;
}

void OutputAlignmentMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  // First get alignment:
  auto aln = block.getAlignment();
  // Format sequence names:
  vector<string> names(aln->getNumberOfSequences());
  for (size_t i = 0; i < aln->getNumberOfSequences(); ++i)
  {
    const MafSequence& mafseq = block.sequence(i);
    if (mafseq.hasCoordinates() && outputCoordinates_)
      names[i] = mafseq.getSpecies() + "-" + mafseq.getChromosome() + "(" + mafseq.getStrand() + ")/" + TextTools::toString(mafseq.start() + 1) + "-" + TextTools::toString(mafseq.stop() + 1);
    else
      names[i] = mafseq.getSpecies();
  }
  aln->setSequenceNames(names, true);
  if (addLDHatHeader_)
    out << aln->getNumberOfSequences() << " " << aln->getNumberOfSites() << " 1" << endl; // We here assume sequences are haploid.
  writer_->writeAlignment(out, *aln);
}
