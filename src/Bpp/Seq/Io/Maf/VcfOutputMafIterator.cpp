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

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

void VcfOutputMafIterator::writeHeader(std::ostream& out) const
{
  time_t t = time(0); // get current time
  struct tm* ct = localtime(&t);
  out << "##fileformat=VCFv4.0" << endl;
  out << "##fileDate=" << (ct->tm_year + 1900) << (ct->tm_mon + 1) << ct->tm_mday << endl;
  out << "##source=Bio++" << endl;  
  out << "##FILTER=<ID=gap,Description=\"At least one sequence contains a gap\">" << endl;
  out << "##FILTER=<ID=unk,Description=\"At least one sequence contains an unresolved character\">" << endl;
  //There are more options in the header that we may want to support...

  //Now write the header line:
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << endl;
}

void VcfOutputMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const
{
  const MafSequence& refSeq = block.getSequenceForSpecies(refSpecies_);
  string chr = refSeq.getChromosome();
  size_t offset = refSeq.start();
  int gap = refSeq.getAlphabet()->getGapCharacterCode();
  string chars = "";
  for (int i = 0; i < static_cast<int>(AlphabetTools::DNA_ALPHABET.getNumberOfTypes()); ++i)
    chars += AlphabetTools::DNA_ALPHABET.intToChar(i);
  VectorSiteContainer sites(block.getAlignment());
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
        filter += ",";
      filter += "unk";
    }
    if (filter == "")
      filter = "PASS";
      
    map<int, size_t> counts;
    SiteTools::getCounts(sites.getSite(i), counts);
    int ref = refSeq[i];
    string alt = "";
    string ac = "";
    for (int x = 0; x < 4; ++x) {
      if (x != ref) {
        size_t f = counts[x];
        if (f > 0) {
          if (alt != "") {
            alt += ",";
            ac += ",";
          }
          alt += TextTools::toString<char>(chars[x]);
          ac += TextTools::toString(f);
        }
      }
    }
    if (ac != "") {
      out << chr << "\t" << (offset + i + 1) << "\t.\t" << chars[refSeq[i]] << "\t" << alt << "\t.\t" << filter << "\tAC=" << ac << endl;
    }
  }
}

