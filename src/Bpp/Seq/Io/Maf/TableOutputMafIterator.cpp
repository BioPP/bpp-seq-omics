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
  // Check for reference species:
  if (block.hasSequenceForSpecies(refSpecies_))
  {
    const auto& refSeq = block.sequenceForSpecies(refSpecies_);
    auto walker = make_unique<SequenceWalker>(refSeq);
    string chr = refSeq.getChromosome();

    // Preprocess data:
    vector<string> seqs;
    for (const string& sp : species_)
    {
      if (block.hasSequenceForSpecies(sp))
      {
        seqs.push_back(block.sequenceForSpecies(sp).toString());
      }
      else
      {
        seqs.push_back("");
      }
    }
    // Loop over all alignment columns:
    for (size_t i = 0; i < block.getNumberOfSites(); ++i)
    {
      string pos = TextTools::toString(walker->getSequencePosition(i));
      *output_ << chr << "\t" << pos;
      for (const string& seq : seqs)
      {
	if (seq != "")
	{
          *output_ << "\t" << seq[i];
	}
	else
	{
          *output_ << "\tNA";
	}
      }
      *output_ << endl;
    }
  } else {
     if (logstream_)
     {
       (*logstream_ << "OUTPUT AS TABLE FILTER: block " << block.getDescription() << " does not contain the reference species.").endLine();
     }
  }
}
