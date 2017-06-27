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

void TableOutputMafIterator::writeBlock_(std::ostream& out, const MafBlock& block)
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

  //first check if there is one (for now we assume that features refer to the chromosome or contig name, with implicit species):
  std::map<std::string, RangeSet<size_t> >::iterator mr = ranges_.find(refSeq.getChromosome());
  if (mr == ranges_.end())
    goto START;
        
  RangeSet<size_t> ranges = mr->second;
  if (completeOnly_)
    ranges.filterWithin(refSeq.getRange(true));
  else  
    ranges.restrictTo(refSeq.getRange(true));
  
  if (!ranges.isEmpty()) {
    //CONTINUE HERE

    SequenceWalker walker(refSeq);
  }
}

