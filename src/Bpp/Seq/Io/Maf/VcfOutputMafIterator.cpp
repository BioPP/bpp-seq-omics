// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "VcfOutputMafIterator.h"

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

void VcfOutputMafIterator::writeHeader_(std::ostream& out) const
{
  time_t t = time(0); // get current time
  struct tm* ct = localtime(&t);
  out << "##fileformat=VCFv4.0" << endl;
  out << "##fileDate=" << (ct->tm_year + 1900) << (ct->tm_mon + 1) << ct->tm_mday << endl;
  out << "##source=Bio++" << endl;
  out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
  out << "##FILTER=<ID=gap,Description=\"At least one sequence contains a gap\">" << endl;
  out << "##FILTER=<ID=unk,Description=\"At least one sequence contains an unresolved character\">" << endl;
  if (genotypes_.size() > 0)
    out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  out << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl;
  // There are more options in the header that we may want to support...

  // Now write the header line:
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
  if (genotypes_.size() > 0)
  {
    out << "\tFORMAT";
    for (size_t i = 0; i < genotypes_.size(); ++i)
    {
      out << "\t" << genotypes_[i][0];
      if (genotypes_[i].size() > 1)
      {
        // Polyploid case
        for (size_t j = 1; j < genotypes_[i].size(); ++j)
        {
          out << "-" << genotypes_[i][j];
        }
      }
    }
  }
  out << endl;
}

void VcfOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block) const
{
  if (block.hasSequenceForSpecies(refSpecies_))
  {
    const MafSequence& refSeq = block.sequenceForSpecies(refSpecies_);
    string chr = refSeq.getChromosome();
    SequenceWalker walker(refSeq);
    size_t offset = refSeq.start();
    int gap = refSeq.alphabet().getGapCharacterCode();
    map<int, string> chars;
    for (int i = (gapAsDeletion_ ? -1 : 0); i < static_cast<int>(AlphabetTools::DNA_ALPHABET->getNumberOfTypes()); ++i)
    {
      chars[i] = AlphabetTools::DNA_ALPHABET->intToChar(i);
    }
    // Where to store genotype information, if any:
    vector<int> gt(genotypes_.size());
    // Now we look all sites for SNPs:
    for (size_t i = 0; i < block.getNumberOfSites(); ++i)
    {
      if (refSeq[i] == gap) // TODO: call indels
        continue;
      string filter = "";
      if (!gapAsDeletion_ && SiteTools::hasGap(block.site(i)))
      {
        filter = "gap";
      }
      if (SymbolListTools::hasUnresolved(block.site(i)))
      {
        if (filter != "")
          filter += ";";
        filter += "unk";
      }
      if (filter == "")
        filter = "PASS";

      map<int, size_t> counts;
      SiteTools::getCounts(block.site(i), counts);
      int ref = refSeq[i];
      string alt = "";
      string ac = "";

      map<int, int> snps;
      int c = 0;
      for (int x = (gapAsDeletion_ ? -1 : 0); x < 4; ++x)
      {
        if (x != ref)
        { 
          size_t f = counts[x];
          if (f > 0)
          {
            if (alt != "")
            {
              alt += ",";
              ac += ",";
            }
            alt += chars[x];
            ac += TextTools::toString(f);
            snps[x] = ++c;
          }
        }
        else
        {
          snps[x] = 0;
        }
      }
      if (ac == "" && outputAll_)
      {
        ac = TextTools::toString(counts[ref]);
      }
      if (ac != "")
      {
        out << chr << "\t" << (offset + walker.getSequencePosition(i) + 1) << "\t.\t" << chars[refSeq[i]] << "\t" << alt << "\t.\t" << filter << "\tAC=" << ac;
        // Write genotpyes:
        if (genotypes_.size() > 0)
        {
          out << "\tGT";
          for (size_t g = 0; g < genotypes_.size(); ++g)
          {
            string geno = "";
            for (auto x: genotypes_[g])
            {
              if (geno != "")
                geno += "|"; // Polyploid
              vector<const MafSequence*> sequences = block.getSequencesForSpecies(x);
              if (sequences.size() == 0)
                geno += (generateDiploids_ ? ".|." : ".");
              else if (sequences.size() > 1)
                throw Exception("VcfOutputMafIterator::writeBlock(). Duplicated sequence for species '" + x + "'.");
              else
              {
                int state = (*sequences[0])[i];
                if (AlphabetTools::DNA_ALPHABET->isUnresolved(state) || (AlphabetTools::DNA_ALPHABET->isGap(state) && !gapAsDeletion_))
                {
                  geno += (generateDiploids_ ? ".|." : ".");
                }
                else
                {
                  geno += TextTools::toString(snps[state]);
                  if (generateDiploids_)
                    geno += "|" + TextTools::toString(snps[state]);
                }
              }
            }
            out << "\t" << geno;
          }
        }
        out << endl;
      }
    }
  } else {
    if (logstream_)
    {
      (*logstream_ << "VCF OUTPUT: block " << block.getDescription() << " does not contain the reference species.").endLine();
    }
  }
}
