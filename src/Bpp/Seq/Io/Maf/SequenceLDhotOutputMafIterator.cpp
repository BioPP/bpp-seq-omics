// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SequenceLDhotOutputMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>

using namespace std;

unique_ptr<MafBlock> SequenceLDhotOutputMafIterator::analyseCurrentBlock_()
{
  auto block = iterator_->nextBlock();
  if (block)
  {
    string chr   = "ChrNA";
    string start = "StartNA";
    string stop  = "StopNA";
    if (block->hasSequenceForSpecies(refSpecies_))
    {
      const MafSequence& refseq = block->sequenceForSpecies(refSpecies_);
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
  return block;
}

void SequenceLDhotOutputMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  // First get alignment:
  unique_ptr<VectorSiteContainer> variableSites = make_unique<VectorSiteContainer>(block.getSequenceNames(), AlphabetTools::DNA_ALPHABET);

  // We first preparse the data:
  // We assume all sequences are distinct:
  size_t nbDistinct = block.getNumberOfSequences();
  size_t nbGenes = block.getNumberOfSequences();
  size_t nbLoci = 0;

  string positions = "";
  for (size_t i = 0; i < block.getNumberOfSites(); ++i)
  {
    const Site& s = block.site(i);
    if (completeOnly_ && !SiteTools::isComplete(s))
    {
      continue;
    }
    unsigned int count = 0;
    int x = -1;
    for (size_t j = 0; j < s.size() && count < 2; ++j)
    {
      if (!AlphabetTools::DNA_ALPHABET->isGap(s[j]) && !AlphabetTools::DNA_ALPHABET->isUnresolved(s[j]))
      {
        if (count == 0)
        {
          // First state found
          count++;
          // We record the state
          x = s[j];
        }
        else
        {
          if (s[j] != x)
          {
            // New state found
            count++;
          }
          // Otherwise, same state as before.
        }
      }
    }
    if (count == 2)
    {
      // At least two alleles (non-gap, non-unresolved) found in this position, so we record it
      positions += " " + TextTools::toString(i + 1);
      nbLoci++;
      auto tmpSite = make_unique<Site>(block.site(i));
      variableSites->addSite(tmpSite);
    }
  }

  // Write header:
  out << "Distinct = " << nbDistinct << endl;
  out << "Genes = " << nbGenes << endl;
  out << "Loci = " << nbLoci << endl;
  out << "K = -4 %4-allele model with Haplotype Alleles specified by A,C,G,T" << endl;

  out << "Positions of loci:" << endl;
  out << positions << endl;

  out << "Haplotypes" << endl;

  for (size_t i = 0; i < block.getNumberOfSequences(); ++i)
  {
    out << variableSites->sequence(i).toString() << " 1" << endl;
  }

  out << "#" << endl;
}
