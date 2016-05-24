//
// File: OutputSequenceLDHotMafIterator.cpp
// Authors: Julien Dutheil
// Created: Thr May 12 2016
//

/*
Copyright or © or Copr. Bio++ Development Team, (2016)

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

#include "SequenceLDhotOutputMafIterator.h"

//From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

//From the STL:
#include <string>
#include <numeric>

using namespace std;

MafBlock* SequenceLDhotOutputMafIterator::analyseCurrentBlock_() throw (Exception)
{
  MafBlock* block = iterator_->nextBlock();
  if (block) {
    string chr   = "ChrNA";
    string start = "StartNA";
    string stop  = "StopNA";
    if (block->hasSequenceForSpecies(refSpecies_)) {
      const MafSequence& refseq = block->getSequenceForSpecies(refSpecies_);
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

void SequenceLDhotOutputMafIterator::writeBlock(std::ostream& out, const MafBlock& block) const {
  //First get alignment:
  const SiteContainer& aln = block.getAlignment();
  unique_ptr<VectorSiteContainer> variableSites(new VectorSiteContainer(aln.getSequencesNames(), &AlphabetTools::DNA_ALPHABET));
  
  //We first preparse the data:
  //We assume all sequences are distinct:
  size_t nbDistinct = aln.getNumberOfSequences();
  size_t nbGenes = aln.getNumberOfSequences();
  size_t nbLoci = 0;

  string positions = "";
  for (size_t i = 0; i < aln.getNumberOfSites(); ++i) {
    const Site& s = aln.getSite(i);
    unsigned int count = 0;
    int x = -1;
    for (size_t j = 0; j < s.size() && count < 2; ++j) {
      if (!AlphabetTools::DNA_ALPHABET.isGap(s[j]) && !AlphabetTools::DNA_ALPHABET.isUnresolved(s[j])) {
        if (count == 0) {
          //First state found
          count++;
          //We record the state
          x = s[j];
        } else {
          if (s[j] != x) {
            //New state found
            count++;
          }
          //Otherwise, same state as before.
        }
      }
    }
    if (count == 2) {
      //At least two alleles (non-gap, non-unresolved) found in this position, so we record it
      positions += " " + TextTools::toString(i + 1);
      nbLoci++;
      variableSites->addSite(aln.getSite(i));
    }
  }

  //Write header:
  out << "Distinct = " << nbDistinct << endl;
  out << "Genes = " << nbGenes << endl;
  out << "Loci = " << nbLoci << endl;
  out << "K = -4 %4-allele model with Haplotype Alleles specified by A,C,G,T" << endl;

  out << "Positions of loci:" << endl;
  out << positions << endl;

  out << "Haplotypes" << endl;

  for (size_t i = 0; i < aln.getNumberOfSequences(); ++i) {
    out << variableSites->getSequence(i).toString() << " 1" << endl;
  }

  out << "#" << endl;
}

