//
// File: TableOutputMafIterator.cpp
// Authors: Julien Dutheil
// Created: Tue Jan 06 2015
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

#include "TableOutputMafIterator.h"

//From bpp-core:
#include <Bpp/Text/TextTools.h>

//From bpp-seq:
#include <Bpp/Seq/SequenceWalker.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>
#include <ctime>

using namespace std;

void TableOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block)
{
  //Check for reference species for coordinates:
  unique_ptr<SequenceWalker> walker;
  bool hasCoordinates = block.hasSequenceForSpecies(refSpecies_);
  string chr = "NA";
  string pos = "NA";
  if (hasCoordinates) {
    const MafSequence& refSeq = block.getSequenceForSpecies(refSpecies_);
    walker.reset(new SequenceWalker(refSeq));
    chr = refSeq.getChromosome();
  }

  //Preprocess data:
  vector<string> seqs;
  for (const string& sp : species_) {
    seqs.push_back(block.getSequenceForSpecies(sp).toString());
  } 
  //Loop over all alignment columns:
  for (size_t i = 0; i < block.getNumberOfSites(); ++i) {
    pos = TextTools::toString(walker->getSequencePosition(i));
    if (hasCoordinates) {
      *output_ << chr << "\t" << pos;
    }
    for (const string& seq : seqs) {
      *output_ << "\t" << seq[i];
    }
    *output_ << endl; 
  }
}

