// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PlinkOutputMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

// From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

void PlinkOutputMafIterator::init_()
{
  nbIndividuals_ = makeDiploids_ ? species_.size() / 2 : species_.size();
  ped_.resize(nbIndividuals_);
  for (size_t i = 0; i < nbIndividuals_; ++i)
  {
    ped_[i] = "FAM001" + colSeparator_ + TextTools::toString(i + 1) + colSeparator_ + "0" + colSeparator_ + "0" + colSeparator_ + "0" + colSeparator_ + TextTools::toString(phenotype_);
  }
}


void PlinkOutputMafIterator::parseBlock_(std::ostream& out, const MafBlock& block)
{
  // Preliminary stuff...
  VectorSiteContainer sites(AlphabetTools::DNA_ALPHABET);
  for (auto& species : species_)
  {
    if (block.hasSequenceForSpecies(species))
    {
      auto tmpSeq = make_unique<Sequence>(block.sequenceForSpecies(species));
      sites.addSequence(tmpSeq->getName(), tmpSeq);
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
  const MafSequence& refSeq = block.sequenceForSpecies(refSpecies_);
  string chr = refSeq.getChromosome();
  if (chr != currentChr_)
  {
    currentChr_ = chr;
    lastPosition_ = 0;
  }
  else
  {
    // Check that block are ordered according to reference sequence:
    if (refSeq.start() < lastPosition_)
      throw Exception("PlinkOutputMafIterator: blocks are not projected according to reference sequence: " + refSeq.getDescription() + "<!>" + TextTools::toString(lastPosition_) + ".");
    lastPosition_ = refSeq.stop();
  }

  string chrStr = chr;
  if (recodeChr_)
  {
    auto i = chrCodes_.find(chr);
    if (i != chrCodes_.end())
    {
      chrStr = TextTools::toString(i->second);
    }
    else
    {
      unsigned int code = currentCode_++;
      chrCodes_[chr] = code;
      chrStr = TextTools::toString(code);
    }
  }

  SequenceWalker walker(refSeq);
  size_t offset = refSeq.start();
  int gap = refSeq.getAlphabet()->getGapCharacterCode();

  // Now we shall scan all sites for SNPs:
  for (size_t i = 0; i < sites.getNumberOfSites(); i++)
  {
    if (refSeq[i] == gap)
      continue;

    // We call SNPs only at position without gap or unresolved characters, and for biallelic sites:
    if (SiteTools::isComplete(sites.site(i)) && SiteTools::getNumberOfDistinctCharacters(sites.site(i)) == 2)
    {
      string pos = "NA";
      if (refSeq[i] != gap)
      {
        pos = TextTools::toString(offset + walker.getSequencePosition(i) + 1);
      }
      string alleles = sites.site(i).toString();
      if (makeDiploids_) {
        for (size_t j = 0; j < nbIndividuals_; ++j)
        {
          ped_[j] += colSeparator_ + TextTools::toString(alleles[2 * j]) + " " + alleles[2 * j + 1];
        }
      } else {
        for (size_t j = 0; j < alleles.size(); ++j)
        {
          ped_[j] += colSeparator_ + TextTools::toString(alleles[j]) + " " + alleles[j];
        }
      }
      // SNP identifier are built as <chr>.<pos>
      out << chrStr << colSeparator_ << chr << "." << pos << colSeparator_;
      if (!map3_)
        out << "0" << colSeparator_; // Add null genetic distance
      out << pos << endl;
    }
  }
}

void PlinkOutputMafIterator::writePedToFile_(ostream& out)
{
  for (auto ped: ped_)
  {
    out << ped << endl;
  }
}
