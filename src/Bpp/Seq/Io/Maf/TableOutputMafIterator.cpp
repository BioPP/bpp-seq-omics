// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "TableOutputMafIterator.h"

// From bpp-core:
#include <Bpp/Text/TextTools.h>

// From bpp-seq:
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

void TableOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block)
{
  // Check for reference species for coordinates:
  unique_ptr<SequenceWalker> walker;
  bool hasCoordinates = block.hasSequenceForSpecies(refSpecies_);
  string chr = "NA";
  string pos = "NA";
  if (hasCoordinates)
  {
    const auto& refSeq = block.sequenceForSpecies(refSpecies_);
    walker.reset(new SequenceWalker(refSeq));
    chr = refSeq.getChromosome();
  }

  // Preprocess data:
  vector<string> seqs;
  for (const string& sp : species_)
  {
    seqs.push_back(block.sequenceForSpecies(sp).toString());
  }
  // Loop over all alignment columns:
  for (size_t i = 0; i < block.getNumberOfSites(); ++i)
  {
    pos = TextTools::toString(walker->getSequencePosition(i));
    if (hasCoordinates)
    {
      *output_ << chr << "\t" << pos;
    }
    for (const string& seq : seqs)
    {
      *output_ << "\t" << seq[i];
    }
    *output_ << endl;
  }
}
