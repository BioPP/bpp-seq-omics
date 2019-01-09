//
// File: VcfOutputMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Jan 05 2013
//

/*
Copyright or © or Copr. Bio++ Development Team, (2010)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "VcfOutputMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Seq/SequenceWithQuality.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

//From the STL:
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
  //There are more options in the header that we may want to support...

  //Now write the header line:
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
  if (genotypes_.size() > 0) {
    out << "\tFORMAT";
    for (size_t i = 0; i < genotypes_.size(); ++i)
      out << "\t" << genotypes_[i];
  }
  out << endl;
}

void VcfOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block) const
{
  const MafSequence& refSeq = block.getSequenceForSpecies(refSpecies_);
  string chr = refSeq.getChromosome();
  SequenceWalker walker(refSeq);
  size_t offset = refSeq.start();
  int gap = refSeq.getAlphabet()->getGapCharacterCode();
  map<int, string> chars;
  for (int i = 0; i < static_cast<int>(AlphabetTools::DNA_ALPHABET.getNumberOfTypes()); ++i)
    chars[i] = AlphabetTools::DNA_ALPHABET.intToChar(i);
  VectorSiteContainer sites(block.getAlignment());
  //Where to store genotype information, if any:
  vector<int> gt(genotypes_.size());
  //Now we look all sites for SNPs:
  for (size_t i = 0; i < sites.getNumberOfSites(); i++) {
    if (refSeq[i] == gap)
      continue;
    string filter = "";
    if (SiteTools::hasGap(sites.getSite(i))) {
      filter = "gap";
    }
    if (SiteTools::hasUnknown(sites.getSite(i))) {
      if (filter != "")
        filter += ";";
      filter += "unk";
    }
    if (filter == "")
      filter = "PASS";
      
    map<int, size_t> counts;
    SiteTools::getCounts(sites.getSite(i), counts);
    int ref = refSeq[i];
    string alt = "";
    string ac = "";
    
    map<int, int> snps;
    int c = 0;
    for (int x = 0; x < 4; ++x) {
      if (x != ref) {
        size_t f = counts[x];
        if (f > 0) {
          if (alt != "") {
            alt += ",";
            ac += ",";
          }
          alt += chars[x];
          ac += TextTools::toString(f);
          snps[x] = ++c;
        }
      } else {
        snps[x] = 0;
      }
    }
    if (ac == "" && outputAll_) {
      ac = TextTools::toString(counts[ref]);
    }
    if (ac != "") {
      out << chr << "\t" << (offset + walker.getSequencePosition(i) + 1) << "\t.\t" << chars[refSeq[i]] << "\t" << alt << "\t.\t" << filter << "\tAC=" << ac;
      //Write genotpyes:
      if (genotypes_.size() > 0) {
        out << "\tGT";
        for (size_t g = 0; g < genotypes_.size(); ++g) {
          vector<const MafSequence*> sequences = block.getSequencesForSpecies(genotypes_[g]);
          if (sequences.size() == 0)
            out << "\t.";
          else if (sequences.size() > 1)
            throw Exception("VcfOutputMafIterator::writeBlock(). Duplicated sequence for species '" + genotypes_[g] + ",.");
          else {
            int state = (*sequences[0])[i];
            if (AlphabetTools::DNA_ALPHABET.isGap(state) || AlphabetTools::DNA_ALPHABET.isUnresolved(state))
              out << (generateDiploids_ ? "\t.|." : "\t.");
            else { 
              out << "\t" << snps[state];
              if (generateDiploids_) 
                out << "|" << snps[state];
            }
          }
        }
      }
      out << endl;
    }
  }
}

