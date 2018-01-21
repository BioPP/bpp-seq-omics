//
// File: PlinkOutputMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Dec 10 2015
//

/*
Copyright or © or Copr. Bio++ Development Team, (2015)

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

#include "PlinkOutputMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/SequenceWithAnnotationTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

void PlinkOutputMafIterator::init_() {
  for (size_t i = 0; i < species_.size(); ++i) {
    ped_[i] = "FAM001\t" + TextTools::toString(i + 1) + "\t0\t0\t0\t0";
  }
}


void PlinkOutputMafIterator::parseBlock_(std::ostream& out, const MafBlock& block)
{
  //Preliminary stuff...
  VectorSiteContainer sites(&AlphabetTools::DNA_ALPHABET);
  for (size_t i = 0; i < species_.size(); ++i) {
    if (block.hasSequenceForSpecies(species_[i])) {
      sites.addSequence(block.getSequenceForSpecies(species_[i]));
      //Note: in case of duplicates, this takes the first sequence.
    } else {
      //Block with missing species are ignored.
      return;
    }
  }
  //Get the reference species for coordinates:
  if (! block.hasSequenceForSpecies(refSpecies_))
    return;
  const MafSequence& refSeq = block.getSequenceForSpecies(refSpecies_);
  string chr = refSeq.getChromosome();
  if (chr != currentChr_) {
    currentChr_ = chr;
    lastPosition_ = 0;
  } else {
    //Check that block are ordered according to reference sequence:
    if (refSeq.start() < lastPosition_)
      throw Exception("PlinkOutputMafIterator: blocks are not projected according to reference sequence: " + refSeq.getDescription() + "<!>" + TextTools::toString(lastPosition_) + ".");
    lastPosition_ = refSeq.stop();
  }

  string chrStr = chr;
  if (recodeChr_) {
    auto i = chrCodes_.find(chr);
    if (i != chrCodes_.end()) {
      chrStr = TextTools::toString(i->second);
    } else {
      unsigned int code = currentCode_++;
      chrCodes_[chr] = code;
      chrStr = TextTools::toString(code);
    }
  }

  SequenceWalker walker(refSeq);
  size_t offset = refSeq.start();
  int gap = refSeq.getAlphabet()->getGapCharacterCode();
  
  //Now we shall scan all sites for SNPs:
  for (size_t i = 0; i < sites.getNumberOfSites(); i++) {
    if (refSeq[i] == gap)
      continue;

    //We call SNPs only at position without gap or unresolved characters, and for biallelic sites:
    if (SiteTools::isComplete(sites.getSite(i)) && SiteTools::getNumberOfDistinctCharacters(sites.getSite(i)) == 2) {
      string pos = "NA";
      if (refSeq[i] != gap) {
        pos = TextTools::toString(offset + walker.getSequencePosition(i) + 1);
      }
      string alleles = sites.getSite(i).toString();
      for (size_t j = 0; j < alleles.size(); ++j) {
        ped_[j] += "\t" + TextTools::toString(alleles[j]) + " " + alleles[j];
      }
      // SNP identifier are built as <chr>.<pos>
      out << chrStr << "\t" << chr << "." << pos << "\t";
      if (!map3_)
        out << "0\t"; //Add null genetic distance
      out << pos << endl;
    }
  }
}

void PlinkOutputMafIterator::writePedToFile_(ostream& out) {
  for (size_t i = 0; i < ped_.size(); ++i)
    out << ped_[i] << endl;
}

