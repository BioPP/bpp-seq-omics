// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MsmcOutputMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

void MsmcOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block)
{
  // Preliminary stuff...

  VectorSiteContainer sites(AlphabetTools::DNA_ALPHABET);
  for (size_t i = 0; i < species_.size(); ++i)
  {
    if (block.hasSequenceForSpecies(species_[i]))
    {
      auto tmpseq = make_unique<Sequence>(block.sequenceForSpecies(species_[i]));
      sites.addSequence(tmpseq->getName(), tmpseq);
      // Note: in case of duplicates, this takes the first sequence.
    }
    else
    {
      // Block with missing species are ignored.
      return;
    }
  }
  // Get the reference species for coordinates:
  if (!block.hasSequenceForSpecies(refSpecies_))
    return;
  const auto& refSeq = block.sequenceForSpecies(refSpecies_);
  string chr = refSeq.getChromosome();
  if (chr != currentChr_)
  {
    currentChr_ = chr;
    nbOfCalledSites_ = 0; // Reset count of called sites.
    lastPosition_ = 0;
  }
  else
  {
    // Check that block are ordered according to reference sequence:
    if (refSeq.start() < lastPosition_)
      throw Exception("MsmcOutputMafIterator: blocks are not projected according to reference sequence: " + refSeq.getDescription() + "<!>" + TextTools::toString(lastPosition_) + ".");
    lastPosition_ = refSeq.stop();
  }

  SequenceWalker walker(refSeq);
  size_t offset = refSeq.start();
  int gap = refSeq.getAlphabet()->getGapCharacterCode();

  // Now we shall scan all sites for SNPs:
  for (size_t i = 0; i < sites.getNumberOfSites(); i++)
  {
    if (refSeq[i] == gap)
      continue;

    // We call SNPs only at position without gap or unresolved characters:
    if (SiteTools::isComplete(sites.site(i)))
    {
      nbOfCalledSites_++;

      if (!SiteTools::isConstant(sites.site(i)))
      {
        string pos = "NA";
        if (refSeq[i] != gap)
        {
          pos = TextTools::toString(offset + walker.getSequencePosition(i) + 1);
        }
        out << chr << "\t" << pos << "\t" << nbOfCalledSites_ << "\t" << sites.site(i).toString() << endl;
        // Reset number of called sites
        nbOfCalledSites_ = 0;
      }
    }
  }
}
