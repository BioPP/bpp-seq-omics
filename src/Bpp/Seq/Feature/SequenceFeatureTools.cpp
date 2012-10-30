//
// File: SequenceFeatureTools.h
// Created by: Julien Dutheil
// Created on: Mon Jul 30 2012
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#include "SequenceFeatureTools.h"

//From bpp-seq:
#include <Bpp/Seq/SequenceTools.h> 
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/AlphabetExceptions.h>

//From STL
#include <vector>

using namespace bpp;

/******************************************************************************/

Sequence* SequenceFeatureTools::extract(const Sequence& seq, const SeqRange& range)
{
  if (range.end() > seq.size())
    throw IndexOutOfBoundsException ("SequenceTools::extract: Invalid upper bound", range.end(), 0, seq.size());
  Sequence* sout = SequenceTools::subseq(seq, range.begin(), range.end() - 1);
  if (range.isNegativeStrand()) {
    SequenceTools::invertComplement(*sout);
  }
  return sout;
}

/******************************************************************************/

unsigned int SequenceFeatureTools::getOrfs(const Sequence& seq, SequenceFeatureSet& featSet) {
  if (! AlphabetTools::isNucleicAlphabet(seq.getAlphabet())) {
    throw AlphabetException("SequenceFeatureTools::getOrfs: Sequence alphabet must be nucleic!", seq.getAlphabet());
  }
  unsigned int orfCpt = 0;
  bpp::StandardCodonAlphabet codonAlpha(dynamic_cast< const NucleicAlphabet* >(seq.getAlphabet()));
  std::vector< std::vector< size_t > > starts(3), stops(3);
  unsigned int phase = 0;
  for (unsigned int p = 0 ; p < seq.size() - 2 ; p++) {
    phase = p % 3;
    if (codonAlpha.isInit(codonAlpha.getCodon(seq.getValue(p), seq.getValue(p + 1), seq.getValue(p + 2)))) {
      starts[phase].push_back(p);
      //std::cerr << "Start: " << p << " (" << phase << ")" << std::endl;
    } else if (codonAlpha.isStop(codonAlpha.getCodon(seq.getValue(p), seq.getValue(p + 1), seq.getValue(p + 2)))) {
      stops[phase].push_back(p);
      //std::cerr << "Stop:  " << p << " (" << phase << ")" << std::endl;
    }
  }
  for (size_t i = 0 ; i < 3 ; ++i) {
    std::vector< size_t >::iterator start(starts[i].begin()), stop(stops[i].begin());
    while (stop != stops[i].end() && start != starts[i].end()) {
      if (*stop < *start) {
        stop++;
      } else {
        orfCpt++;
        //std::cerr << "ORF:  " << *start << " - " << *stop + 2 << " (" << i << ")" << std::endl;
        bpp::BasicSequenceFeature feat("", seq.getName(), "Bio++", "CDS", *start, *stop + 2, '+');
        featSet.addFeature(feat);
        start++;
      }
    }
  }
  return orfCpt;
}

/******************************************************************************/
